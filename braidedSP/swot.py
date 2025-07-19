# Typing imports
from dataclasses import dataclass
import os
from tqdm import tqdm
from datetime import datetime
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import rasterio
from rasterio.transform import rowcol


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


    def load_pixc(self, extraction_mask, mask_transform, mask_crs, centerline_mask=[], engine='h5netcdf'):

        # read in initial pixc data, clip if there is an initial buffered centerline mask provided
        pixc_gdf = _readPIXC(self.path, gdf_buffered=centerline_mask, engine=engine)

        # Trim data to simplified extraction mask
        self.gdf = _trim2mask_general(pixc_gdf, extraction_mask, mask_transform, mask_crs) # Trimming to largest component mask




# ------------------------------------------------------------------------------------------
# helper functions
def _readPIXC(filename, gdf_buffered, engine='h5netcdf'):

    # If no centerline buffer to trim to, set gdf_buffered = []
    nc = xr.open_mfdataset(filename, group = 'pixel_cloud', engine=engine)
    
    # Set crs of nc file
    nc = nc.rio.write_crs("EPSG:4326", inplace=True)

    # Set the nc to a geopandas dataset
    class_flat = nc.classification.values.ravel()
    lon_flat = nc.longitude.values.ravel()[class_flat == 4]
    lat_flat = nc.latitude.values.ravel()[class_flat == 4]
    height = nc.height.values.ravel()[class_flat == 4]
    geoid = nc.geoid.values.ravel()[class_flat == 4]
    water_frac = nc.water_frac.values.ravel()[class_flat == 4]
    phase_noise_std = nc.phase_noise_std.values.ravel()[class_flat == 4]
    dheight_dphase = nc.dheight_dphase.values.ravel()[class_flat == 4]
    sig0 = nc.sig0.values.ravel()[class_flat == 4]
    heightEGM = height-geoid

    # create geodataframe
    gdf = gpd.GeoDataFrame(pd.DataFrame({"height":height,"heightEGM":heightEGM,"geoid":geoid,"lat":lat_flat,"lon":lon_flat,"class":class_flat[class_flat == 4],"water_frac":water_frac,"phase_noise_std":phase_noise_std,"dheight_dphase":dheight_dphase,"sig0":sig0}),geometry=gpd.points_from_xy(lon_flat,lat_flat))
    gdf.set_crs(epsg=4326, inplace=True)

    ## Trim data to buffer around imported centerline
    if gdf_buffered:
        gdf = gpd.clip(gdf, gdf_buffered)

    return gdf

def _trim2mask_general(pixc_gdf, extraction_mask, mask_transform, mask_crs):
    # - Function to trim pixel cloud data to a binary mask (numpy array of 1s and 0s)
    # - Water mask generated from S2 optical imagery and imported as a geotiff

    # Open the GeoTIFF file
    # with rasterio.open(mask_tiff_filename) as water_mask:
    #     mask_transform = water_mask.transform  # Get the affine transform
    #     mask_crs = water_mask.crs  # Get the CRS (coordinate reference system)

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
        if 0 <= row < extraction_mask.shape[0] and 0 <= col < mask_transform.shape[1]:
            sampled_values.append(mask_transform[row, col])
        else:
            sampled_values.append(0) 

    pixc_gdf['water_masked'] = np.array(sampled_values) == 1

    # Filter the GeoDataFrame to keep only points within the mask
    trimmed_gdf = pixc_gdf[pixc_gdf['water_masked']].reset_index(drop=True)
    
    return trimmed_gdf

"""
def _get_RiverSP(trimmed_pixc_gdf, cl_gen, hemi, season, pixcdate, ileID, figdir, odir, buffer_conv=50):

    # Function for extracting a 'riverSP' product over all branches of the river centerline
    # For rotating the buffer
    theta = np.radians(90)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c,-s), (s, c))) #2D rotation matrix

    riverSP_gdf = gpd.GeoDataFrame() # initialize dataframe
    iter_max = 5

    utm_crs = cl_gen.estimate_utm_crs()

    for channelID in cl_gen.branch_id:

        swot_break = False

        # Select channel
        sel_channel = cl_gen.loc[cl_gen.branch_id == channelID]

        # convert to UTM
        sel_channel_utm = sel_channel.to_crs(utm_crs)
        coords_list = list(sel_channel_utm.geometry.iloc[0].coords)
        x_coords = [x for x,y in coords_list]
        y_coords = [y for x,y in coords_list]

        # extract point locations every 100 m along the line
        sel_channel_utm["point_list"] = sel_channel_utm.apply(lambda x: create_points(row=x, point_separation=200), axis=1)

        #Create a point dataframe by exploding the point list into individual points/rows
        points_df = sel_channel_utm.explode(column="point_list")
        points_df = points_df.set_geometry("point_list").drop(columns="geometry").rename_geometry("geometry").set_crs(sel_channel_utm.crs)
        points_df = points_df.reset_index(drop=True)

        if len(points_df) <= 1:
            print('Channel contains only one point....skipping channel No. ',channelID)
            continue

        # Extract data from SWOT PIXC based on large manual buffer (larger than any expected channel width)
        # Buffer width and transect length dependant on maximum expected channel width (knowledge from visual inspection of S2 masks) for each season
        if season == 'LF':
            init_buffer_width = 2000
            transect_len = 2000
        else:
            init_buffer_width = 3000
            transect_len = 3000

        buffer_width_prev = init_buffer_width*100 # set large previous buffer width to begin with
        buffer_width = init_buffer_width

        iter = 1
        while abs(buffer_width - buffer_width_prev) > buffer_conv and swot_break == False and iter <= iter_max: # condition for buffer width convergence

            #print('Updating buffer width:',buffer_width)

            #print('Extracting SWOT data within channel buffer...')
            tile_figdir = figdir+'/'+tileID+'/generated_nodes/'
            sub_swot = trim_to_one_channel(sel_channel, trimmed_pixc_gdf,hemi,tile_figdir,pixcdate,channelID,buffer_width=buffer_width,savePlot=True) # Extracts SWOT data for one channel
            sub_swot_utm = sub_swot.to_crs(utm_crs)

            # if sub_swot contains no data...break and continue
            if len(sub_swot.heightEGM) <= 1:
                swot_break = True
                continue

            # Give each PIXC cloud point a point geometry
            mpt = MultiPoint([shape(row['geometry']) for _, row in sub_swot_utm.iterrows()])

            ratioval = 0.01
            if float(sel_channel_utm.geometry.length) > 3000: # change ratio value for SWOT extent polygon determination based on selected channel reach length
                ratioval = 0.05

            shapeout = shapely.concave_hull(mpt, ratio=ratioval,allow_holes=True)
            swot_boundary = gpd.GeoDataFrame(geometry=[shapeout],crs=utm_crs)

            # save transects to a gdf for plotting
            transects_gdf = gpd.GeoDataFrame() # initialize dataframe
            intersects_gdf = gpd.GeoDataFrame() # initialize dataframe
            extract_polys = gpd.GeoDataFrame() # initialize dataframe
            riverSP_gdf_oneChannel = gpd.GeoDataFrame() # initialize dataframe for saving single channels


            #print('Extracting transects, widths and heights...')
            for idx in points_df.index:
                coord = list(points_df.geometry.iloc[idx].coords)
                if idx == 0:
                    coord_b4 = coord
                    coord_af = list(points_df.geometry.iloc[idx+1].coords)

                if idx == len(points_df.index) - 1:
                    coord_b4 = list(points_df.geometry.iloc[idx-1].coords)
                    coord_af = coord


                if idx != len(points_df.index) - 1 and idx != 0:
                    coord_b4 = list(points_df.geometry.iloc[idx-1].coords)
                    coord_af = list(points_df.geometry.iloc[idx+1].coords)

                coord_b4 = [coord_b4[0][0],coord_b4[0][1]]
                coord_af = [coord_af[0][0],coord_af[0][1]]

                # Compute vector and rotated vector
                vec = calculate_vector_general(coord_b4,coord_af)
                rot_vec = np.dot(vec, R)

                pt1 = np.array(coord) - transect_len*rot_vec
                pt2 = np.array(coord) + transect_len*rot_vec
                line = LineString([Point(pt1),Point(pt2)])

                # plot transect
                transect_temp = gpd.GeoDataFrame(geometry=[line],crs=points_df.crs)
                transects_gdf = pd.concat([transects_gdf,transect_temp])

                # Find intersection line with swot polygon
                inter_geom = shapely.intersection(swot_boundary,transect_temp)
                intersects_gdf = pd.concat([intersects_gdf,inter_geom])
                
                # Extract 100 m buffer poly around line (for med H estimation)
                extract_poly = inter_geom.geometry.iloc[0].buffer(100,cap_style="flat")
                extract_poly_gdf = gpd.GeoDataFrame(geometry=[extract_poly],crs=utm_crs)
                extract_polys = pd.concat([extract_polys,extract_poly_gdf])

                sub_swot_utm['inPoly']=extract_poly_gdf.geometry.iloc[0].contains(sub_swot_utm.geometry)
                subset = sub_swot_utm[sub_swot_utm['inPoly']==True]

                lineW =inter_geom.geometry.length

                # Compute median if data exists within subset
                if len(subset) > 0:
                    medH = np.nanmedian(subset.heightEGM)
                    medH_ell = np.nanmedian(subset.height)
                    medGeoid = np.nanmedian(subset.geoid)
                    med_water_frac = np.nanmedian(subset.water_frac)
                    med_phstd = np.nanmedian(subset.phase_noise_std)
                    med_dhdp = np.nanmedian(subset.dheight_dphase)
                    med_sig0 = np.nanmedian(subset.sig0) # changed these to nanmedian
                else:
                    medH = np.nan
                    medH_ell = np.nan
                    medGeoid = np.nan
                    med_water_frac = np.nan
                    med_phstd = np.nan
                    med_dhdp = np.nan
                    med_sig0 = np.nan

                # get line width 

                
                # print('median H:',medH)
                # print('width:',lineW)

                # Save to RiverSP GDF with point geometry
                riverSP_temp = gpd.GeoDataFrame(geometry=[Point(coord)],crs=utm_crs)
                riverSP_temp['heightEGM_med'] = medH
                riverSP_temp['width'] = lineW
                riverSP_temp['channelID'] = channelID
                nodeID = str(channelID) + str(idx)
                riverSP_temp['nodeID'] = nodeID
                riverSP_temp['heightEll_med'] = medH_ell
                riverSP_temp['geoid_med'] = medGeoid
                riverSP_temp['wf_med'] = med_water_frac
                riverSP_temp['phstd_med'] = med_phstd
                riverSP_temp['dhdp_med'] = med_dhdp
                riverSP_temp['sig0_med'] = med_sig0

                riverSP_gdf_oneChannel = pd.concat([riverSP_gdf_oneChannel,riverSP_temp])
            
            iter = iter + 1

            # Update buffer width condition
            medW = np.median(riverSP_gdf_oneChannel.width)
            buffer_width_prev = buffer_width
            buffer_width = np.round(medW*1.1) # buffer with 10% extra of the median width
            #print('Median width of channel:',medW)


        if iter == iter_max:
            print('Buffer width not converged, max iteration reached.')
            

        if swot_break == True:
            print('No SWOT data around channel...skipping channel: ',channelID)
            continue
        
        #print('Final channel width:',medW)
        #print('Final buffer width:',buffer_width_prev)
        riverSP_gdf = pd.concat([riverSP_gdf,riverSP_gdf_oneChannel])


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


        isExist = os.path.exists(figdir+'/'+tileID+'/generated_nodes/')
        if not isExist:
            os.makedirs(figdir+'/'+tileID+'/generated_nodes/')
        plt.savefig(figdir+'/'+tileID+'/generated_nodes/'+str(pixcdate)+'_genNodes_ch'+str(channelID)+'.png')
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
        plt.savefig(figdir+'/'+tileID+'/generated_nodes/'+str(pixcdate)+'_SWOTBOUNDARY_ch'+str(channelID)+'.png')
        plt.close()

        # print('RiverSP for channel No. ',channelID)
        # print(riverSP_gdf.loc[riverSP_gdf.channelID == channelID])

        # Save the SWOT boundary (buffer) used here to extract PIXC data
        testfile = glob.glob(odir+'riverSP_out/'+tileID+'/'+pixcdate+'_'+str(channelID)+'_subswot.geojson')
        if not testfile:    
            print('Saving boundary of selected SWOT data...')
            sub_swot.to_file(odir+'riverSP_out/'+tileID+'/'+pixcdate+'_'+str(channelID)+'_subswot.geojson')
            print('SUCCESS!')


    return riverSP_gdf
"""