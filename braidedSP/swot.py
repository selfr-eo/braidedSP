# Typing imports
from dataclasses import dataclass
import os
from datetime import datetime
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import MultiPoint, LineString, Point, shape
import shapely
from osgeo import gdal  # import before rasterio
from rasterio.transform import rowcol
import rioxarray
import xarray as xr
from tqdm import tqdm

from braidedSP import tools


@dataclass
class SWOT:
    river_name: str
    date: datetime
    path: str
    odir: str

    def __post_init__(self):
        # set empty attributes to be filled
        self.gdf = None

        # make sure file exists
        if not os.path.exists(self.path):
            raise FileNotFoundError(
                "A mask can only be created on an existing raster file."
            )

        # set up export dirs
        if not os.path.exists(self.odir):
            os.makedirs(self.odir)

        # set up fig directory
        self.figdir = os.path.join(self.odir, "figs")
        if not os.path.exists(self.figdir):
            os.makedirs(self.figdir)

    def load_pixc(
        self,
        extraction_mask,
        mask_transform,
        mask_crs,
        centerline_mask=[],
        classes=['open_water'],
        engine="h5netcdf",
    ):
        # read in initial pixc data, clip if there is an initial buffered centerline mask provided
        pixc_gdf = _readPIXC(self.path, gdf_buffered=centerline_mask, classes=classes, engine=engine)

        # Trim data to simplified extraction mask
        self.gdf = _trim2mask_general(
            pixc_gdf, extraction_mask, mask_transform, mask_crs
        )  # Trimming to largest component mask

        # always save timmed data in local crs
        local_crs = self.gdf.estimate_utm_crs()
        self.gdf = self.gdf.to_crs(local_crs)

    def cal_SP(self, cl_gen):
        self.sp = _get_RiverSP(
            self.gdf, cl_gen, buffer_conv=50, init_buffer_width=2000, transect_len=2000
        )


# ------------------------------------------------------------------------------------------
# helper functions
def _translateClass(flag_meanings):

    class_dict = {
        'land':1,
        'land_near_water':2,
        'water_near_land':3,
        'open_water':4,
        'dark_water':5,
        'low_coh_water_near_land':6,
        'open_low_coh_water':7
    }
    flag_values = [class_dict[flag_meaning] for flag_meaning in flag_meanings]

    return flag_values

def _readPIXC(filename, gdf_buffered, classes=['open_water'], engine="h5netcdf"):
    # If no centerline buffer to trim to, set gdf_buffered = []
    nc = xr.open_mfdataset(filename, group="pixel_cloud", engine=engine)

    # Set crs of nc file
    nc = nc.rio.write_crs("EPSG:4326", inplace=True)

    # select based on desired classification
    class_flat = nc.classification.values.ravel()
    flag_values = _translateClass(classes)
    class_condition = np.isin(class_flat, flag_values)

    # Set the nc to a geopandas dataset
    lon_flat = nc.longitude.values.ravel()[class_condition]
    lat_flat = nc.latitude.values.ravel()[class_condition]
    height = nc.height.values.ravel()[class_condition]
    geoid = nc.geoid.values.ravel()[class_condition]
    solid_earth_tide = nc.solid_earth_tide.values.ravel()[class_condition]
    load_tide = nc.load_tide_fes.values.ravel()[class_condition]
    pole_tide = nc.pole_tide.values.ravel()[class_condition]
    water_frac = nc.water_frac.values.ravel()[class_condition]
    phase_noise_std = nc.phase_noise_std.values.ravel()[class_condition]
    dheight_dphase = nc.dheight_dphase.values.ravel()[class_condition]
    sig0 = nc.sig0.values.ravel()[class_condition]

    # correction for solid earth/load/pole tide effects (e.g., see SWOT User Handbook, section 3.1.25)
    heightEGM = height - geoid - solid_earth_tide - load_tide - pole_tide

    # create geodataframe
    data = {
        "height": height,
        "heightEGM": heightEGM,
        "geoid": geoid,
        "lat": lat_flat,
        "lon": lon_flat,
        "class": class_flat[class_condition],
        "water_frac": water_frac,
        "phase_noise_std": phase_noise_std,
        "dheight_dphase": dheight_dphase,
        "sig0": sig0,
    }
    gdf = gpd.GeoDataFrame(
        pd.DataFrame(data), geometry=gpd.points_from_xy(lon_flat, lat_flat)
    )
    gdf.set_crs(epsg=4326, inplace=True)

    ## Trim data to buffer around imported centerline
    if gdf_buffered:
        gdf = gpd.clip(gdf, gdf_buffered)

    return gdf


def _trim2mask_general(pixc_gdf, extraction_mask, mask_transform, mask_crs):
    # - Function to trim pixel cloud data to a binary mask (numpy array of 1s and 0s)
    # - Water mask generated from S2 optical imagery and imported as a geotiff

    # Reproject the GeoDataFrame to match the mask's CRS
    if pixc_gdf.crs != mask_crs:
        pixc_gdf = pixc_gdf.to_crs(mask_crs)

    # Extract point coordinates
    point_coords = [(geom.x, geom.y) for geom in pixc_gdf.geometry]

    # Convert coordinates to array indices
    indices = [rowcol(mask_transform, x, y) for x, y in point_coords]

    # Check for valid indices and sample mask values
    sampled_values = []
    for row, col in indices:
        if 0 <= row < extraction_mask.shape[0] and 0 <= col < extraction_mask.shape[1]:
            sampled_values.append(extraction_mask[row, col])
        else:
            sampled_values.append(0)

    pixc_gdf["water_masked"] = np.array(sampled_values) == 1

    # Filter the GeoDataFrame to keep only points within the mask
    trimmed_gdf = pixc_gdf[pixc_gdf["water_masked"]].reset_index(drop=True)

    return trimmed_gdf


def _get_RiverSP(
    trimmed_pixc_gdf, cl_gen, buffer_conv=50, init_buffer_width=2000, transect_len=2000
):
    # Function for extracting a 'riverSP' product over all branches of the river centerline
    # For rotating the buffer
    theta = np.radians(90)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c, -s), (s, c)))  # 2D rotation matrix

    riverSP_gdf = gpd.GeoDataFrame()  # initialize dataframe
    iter_max = 5

    utm_crs = cl_gen.estimate_utm_crs()

    for channelID in tqdm(cl_gen.branch_id, desc="Processing riverSP for channels"):
        swot_break = False

        # Select channel
        sel_channel = cl_gen.loc[cl_gen.branch_id == channelID].reset_index(drop=True)
        sel_channel_utm = sel_channel.to_crs(utm_crs)

        # extract point locations every 100 m along the line
        sel_channel_utm["point_list"] = sel_channel_utm.apply(
            lambda x: tools.create_points(row=x, point_separation=200), axis=1
        )

        # Create a point dataframe by exploding the point list into individual points/rows
        points_df = sel_channel_utm.explode(column="point_list")
        points_df = (
            points_df.set_geometry("point_list")
            .drop(columns="geometry")
            .rename_geometry("geometry")
            .set_crs(sel_channel_utm.crs)
        )
        points_df = points_df.reset_index(drop=True)
        if len(points_df) <= 1:
            # print('Channel contains only one point....skipping channel No. ', channelID)
            continue

        # Extract data from SWOT PIXC based on large manual buffer (larger than any expected channel width)
        buffer_width_prev = (
            init_buffer_width * 100
        )  # set large previous buffer width to begin with
        buffer_width = init_buffer_width

        iter = 1
        while (
            abs(buffer_width - buffer_width_prev) > buffer_conv
            and not swot_break
            and iter <= iter_max
        ):  # condition for buffer width convergence
            sub_swot = tools.trim_to_one_channel(
                sel_channel, trimmed_pixc_gdf, buffer_width=buffer_width
            )  # Extracts SWOT data for one channel
            sub_swot_utm = sub_swot.to_crs(utm_crs)

            # if sub_swot contains no data...break and continue
            if len(sub_swot) <= 1:
                swot_break = True
                continue

            # print('processing')

            # Give each PIXC cloud point a point geometry
            mpt = MultiPoint(
                [shape(row["geometry"]) for _, row in sub_swot_utm.iterrows()]
            )

            ratioval = 0.01
            if (
                float(sel_channel_utm.loc[0, "geometry"].length) > 3000
            ):  # change ratio value for SWOT extent polygon determination based on selected channel reach length
                ratioval = 0.05

            shapeout = shapely.concave_hull(mpt, ratio=ratioval, allow_holes=True)
            swot_boundary = gpd.GeoDataFrame(geometry=[shapeout], crs=utm_crs)

            # save transects to a gdf for plotting
            transects_gdf = gpd.GeoDataFrame()  # initialize dataframe
            intersects_gdf = gpd.GeoDataFrame()  # initialize dataframe
            extract_polys = gpd.GeoDataFrame()  # initialize dataframe
            riverSP_gdf_oneChannel = (
                gpd.GeoDataFrame()
            )  # initialize dataframe for saving single channels

            # print('Extracting transects, widths and heights...')
            for idx in points_df.index:
                coord = list(points_df.geometry.iloc[idx].coords)
                if idx == 0:
                    coord_b4 = coord
                    coord_af = list(points_df.geometry.iloc[idx + 1].coords)

                if idx == len(points_df.index) - 1:
                    coord_b4 = list(points_df.geometry.iloc[idx - 1].coords)
                    coord_af = coord

                if idx != len(points_df.index) - 1 and idx != 0:
                    coord_b4 = list(points_df.geometry.iloc[idx - 1].coords)
                    coord_af = list(points_df.geometry.iloc[idx + 1].coords)

                coord_b4 = [coord_b4[0][0], coord_b4[0][1]]
                coord_af = [coord_af[0][0], coord_af[0][1]]

                # Compute vector and rotated vector
                vec = tools.calculate_vector_general(coord_b4, coord_af)
                rot_vec = np.dot(vec, R)

                pt1 = np.array(coord) - transect_len * rot_vec
                pt2 = np.array(coord) + transect_len * rot_vec
                line = LineString([Point(pt1), Point(pt2)])

                # plot transect
                transect_temp = gpd.GeoDataFrame(geometry=[line], crs=points_df.crs)
                transects_gdf = pd.concat([transects_gdf, transect_temp])

                # Find intersection line with swot polygon
                inter_geom = shapely.intersection(swot_boundary, transect_temp)
                intersects_gdf = pd.concat([intersects_gdf, inter_geom])

                # Extract 100 m buffer poly around line (for med H estimation)
                extract_poly = inter_geom.geometry.iloc[0].buffer(100, cap_style="flat")
                extract_poly_gdf = gpd.GeoDataFrame(
                    geometry=[extract_poly], crs=utm_crs
                )
                extract_polys = pd.concat([extract_polys, extract_poly_gdf])

                sub_swot_utm["inPoly"] = extract_poly_gdf.geometry.iloc[0].contains(
                    sub_swot_utm.geometry
                )
                subset = sub_swot_utm[sub_swot_utm["inPoly"].values]

                lineW = inter_geom.geometry.length

                # Compute median if data exists within subset
                if len(subset) > 0:
                    medH = np.nanmedian(subset.heightEGM)
                    medH_ell = np.nanmedian(subset.height)
                    medGeoid = np.nanmedian(subset.geoid)
                    med_water_frac = np.nanmedian(subset.water_frac)
                    med_phstd = np.nanmedian(subset.phase_noise_std)
                    med_dhdp = np.nanmedian(subset.dheight_dphase)
                    med_sig0 = np.nanmedian(subset.sig0)  # changed these to nanmedian
                else:
                    medH = np.nan
                    medH_ell = np.nan
                    medGeoid = np.nan
                    med_water_frac = np.nan
                    med_phstd = np.nan
                    med_dhdp = np.nan
                    med_sig0 = np.nan

                # Save to RiverSP GDF with point geometry
                riverSP_temp = gpd.GeoDataFrame(geometry=[Point(coord)], crs=utm_crs)
                riverSP_temp["heightEGM_med"] = medH
                riverSP_temp["width"] = lineW
                riverSP_temp["channelID"] = channelID
                nodeID = str(channelID) + str(idx)
                riverSP_temp["nodeID"] = nodeID
                riverSP_temp["heightEll_med"] = medH_ell
                riverSP_temp["geoid_med"] = medGeoid
                riverSP_temp["wf_med"] = med_water_frac
                riverSP_temp["phstd_med"] = med_phstd
                riverSP_temp["dhdp_med"] = med_dhdp
                riverSP_temp["sig0_med"] = med_sig0

                # drop rows with nan valeus
                riverSP_temp = riverSP_temp.dropna(how='all', axis=0)
                riverSP_temp = riverSP_temp.dropna(how='all', axis=1)

                if len(riverSP_temp) > 0:
                    riverSP_gdf_oneChannel = pd.concat(
                        [riverSP_gdf_oneChannel, riverSP_temp]
                    )

            iter = iter + 1

            # Update buffer width condition
            medW = np.median(riverSP_gdf_oneChannel.width)
            buffer_width_prev = buffer_width
            buffer_width = np.round(
                medW * 1.1
            )  # buffer with 10% extra of the median width
            # print('Median width of channel:',medW)

        if iter == iter_max:
            # print('Buffer width not converged, max iteration reached.')
            pass

        if swot_break:
            # print('No SWOT data around channel...skipping channel: ',channelID)
            continue

        # print('Final channel width:',medW)
        # print('Final buffer width:',buffer_width_prev)
        riverSP_temp = riverSP_temp.dropna(how='all', axis=0)
        riverSP_temp = riverSP_temp.dropna(how='all', axis=1)

        if len(riverSP_gdf_oneChannel) > 0:
            riverSP_gdf = pd.concat([riverSP_gdf, riverSP_gdf_oneChannel])

        """
        # Plot selected channel centerline, extracted points, and transects for each channel
        ax = sel_channel_utm.plot(figsize=(5,5))
        points_df.plot(ax=ax, zorder=2, color="black", markersize=50, marker=".")
        transects_gdf.plot(ax=ax, zorder=2, color="red")
        intersects_gdf.plot(ax=ax, zorder=2, color="yellow")
        sub_swot_utm.plot(ax=ax, column='heightEGM', markersize=0.1, legend=True,cmap='Spectral')
        extract_polys.plot(ax=ax, zorder=2,color='blue',alpha=0.4)
        ax.set_title('Ch. '+str(channelID))
        # Remove x and y tick marks
        ax.set_xticks([])
        ax.set_yticks([])


        isExist = os.path.exists(figdir+'/generated_nodes/')
        if not isExist:
            os.makedirs(figdir+'/generated_nodes/')
        plt.savefig(figdir+'/generated_nodes/'+str(pixcdate)+'_genNodes_ch'+str(channelID)+'.png')
        plt.close()

        # Save also the swot generated polygon
        ax = sel_channel_utm.plot(figsize=(5,5))
        swot_boundary.plot(ax=ax)
        points_df.plot(ax=ax, zorder=2, color="black", markersize=50, marker=".")
        transects_gdf.plot(ax=ax, zorder=2, color="red")
        intersects_gdf.plot(ax=ax, zorder=2, color="yellow")
        ax.set_title('Ch. '+str(channelID))
        # Remove x and y tick marks
        ax.set_xticks([])
        ax.set_yticks([])
        plt.savefig(figdir+'/generated_nodes/'+str(pixcdate)+'_SWOTBOUNDARY_ch'+str(channelID)+'.png')
        plt.close()

        # print('RiverSP for channel No. ',channelID)
        # print(riverSP_gdf.loc[riverSP_gdf.channelID == channelID])

        # Save the SWOT boundary (buffer) used here to extract PIXC data
        
        testfile = glob.glob(odir+'riverSP_out/'+pixcdate+'_'+str(channelID)+'_subswot.geojson')
        if not testfile:    
            print('Saving boundary of selected SWOT data...')
            sub_swot.to_file(odir+'riverSP_out/'+tileID+'/'+pixcdate+'_'+str(channelID)+'_subswot.geojson')
            print('SUCCESS!')
        """

    return riverSP_gdf
