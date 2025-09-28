import geopandas as gpd
import pandas as pd
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import sys

from shapely.geometry import LineString, Point
from shapely import buffer
from skimage.morphology import label, binary_dilation

# RiverSP processing
from shapely.geometry import MultiPoint, shape
import shapely

sys.path.append("../tools")
import procBraided as pc
import skeletonize_func as skel

import warnings
warnings.filterwarnings("ignore")

# Transient functions

def calculate_vector_general(coord1, coord2):
    # Calculate the direction vector for a LineString at the given endpoint.
    vec = np.array(coord2) - np.array(coord1)
    return vec / np.linalg.norm(vec)  # Normalize

# Function from stack exchange: https://gis.stackexchange.com/questions/177583/interpolating-every-x-distance-along-line-in-shapely
def create_points(row, point_separation):
    # Can be used for many different linestrings!
    geom = row.geometry
    #For each geometry, create a point along it for each distance in the range from 0 to line length, with an interval
    point_list = [geom.interpolate(distance=x) for x in np.arange(start=0, stop=geom.length, step=point_separation)]
    return point_list



def trim_to_one_channel(sel_cl_gdf, swot_gdf,hemi,tile_figdir,pixcdate,channelID,buffer_width=1000,savePlot=False):

    # # Buffer around channel (= 1000 unless known riv width passed in)
    # buffer_geom = buffer(sel_cl_gdf.geometry.unary_union, buffer_width)

    # # Trim SWOT PIXC to buffer
    # clipped_gdf = gpd.clip(swot_gdf, buffer_geom) # Use 'trimmed_gdf' to use water mask trimmed data
    # clipped_gdf = clipped_gdf.copy()

    # Buffer around channel (= 1000 unless known riv width passed in)
    gdf_buffered = pc.get_width_dependent_buffer(sel_cl_gdf,'EPSG:4326',hemi,buffer_width,cap_val='flat')

    # Trim SWOT PIXC to buffer
    clipped_gdf = gpd.clip(swot_gdf, gdf_buffered) # Use 'trimmed_gdf' to use water mask trimmed data
    clipped_gdf = clipped_gdf.copy()

    if len(clipped_gdf) > 0:

        # project SWOT PIXC to channel centerline
        clipped_gdf['dist'] = pc.projectToCenterline(sel_cl_gdf,clipped_gdf[['geometry']],hemi)
        clipped_gdf['channelID'] = int(sel_cl_gdf.branch_id)

    # if savePlot==True and len(clipped_gdf) > 0:
    #     fig, ax = plt.subplots(figsize=(10, 10))
    #     swot_gdf.plot(ax=ax, color='purple',markersize=0.1)
    #     clipped_gdf.plot(ax=ax, color='yellow',markersize=0.1)
    #     sel_cl_gdf.plot(ax=ax, color='blue', label='Centerline', linewidth=0.5)
    #     gdf_buffered.plot(ax=ax, color='lightblue', alpha=0.5, label='Buffered Area')
    #     ax.set_title('Channel'+str(channelID)+', Buffer width:'+str(buffer_width),fontsize=15)
    #     # Remove x and y tick marks
    #     ax.set_xticks([])
    #     ax.set_yticks([])
    #     plt.savefig(tile_figdir+str(pixcdate)+'_buffer_ch'+str(channelID)+'.png')
    #     #plt.savefig(figdir+'/'+tileID+'/generated_nodes/'+str(pixcdate)+'_genNodes_ch'+str(channelID)+'.png')
    #     plt.close()

    return clipped_gdf



def get_RiverSP(trimmed_pixc_gdf,cl_gen,hemi,season,pixcdate,tileID,figdir,odir,buffer_conv=50):

    # Function for extracting a 'riverSP' product over all branches of the river centerline
    # For rotating the buffer
    theta = np.radians(90)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c,-s), (s, c))) #2D rotation matrix

    riverSP_gdf = gpd.GeoDataFrame() # initialize dataframe
    iter_max = 5

    for channelID in cl_gen.branch_id:
    #for channelID in [21,26]:
        print('Running channel...',channelID)

        swot_break = False

        #print('.......running width extraction for ch. ',channelID)
        # Select channel
        sel_channel = cl_gen.loc[cl_gen.branch_id == channelID]

        # convert to UTM
        centroid = sel_channel.geometry.unary_union.centroid
        utm_zone = int((centroid.x + 180) // 6) + 1
        if hemi == 'north':
            utm_crs = f'EPSG:{32600 + utm_zone}'
        elif hemi == 'south':
            utm_crs = f'EPSG:{32700 + utm_zone}'
        else:
            raise ValueError("Invalid hemisphere. Use 'north' or 'south'.")

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

            shapeout = shapely.concave_hull(mpt, ratio=ratioval, allow_holes=True)
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
        plt.savefig(figdir+'/'+tileID+'/generated_nodes/'+str(pixcdate)+'_genNodes_ch'+str(channelID)+'_ow_wnl.png')
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
        plt.savefig(figdir+'/'+tileID+'/generated_nodes/'+str(pixcdate)+'_SWOTBOUNDARY_ch'+str(channelID)+'_ow_wnl.png')
        plt.close()

        # print('RiverSP for channel No. ',channelID)
        # print(riverSP_gdf.loc[riverSP_gdf.channelID == channelID])

        # Save the SWOT boundary (buffer) used here to extract PIXC data
        testfile = glob.glob(odir+'riverSP_out/'+tileID+'/'+pixcdate+'_'+str(channelID)+'_subswot_ow_wnl.geojson')
        if not testfile:    
            print('Saving boundary of selected SWOT data...')
            sub_swot.to_file(odir+'riverSP_out/'+tileID+'/'+pixcdate+'_'+str(channelID)+'_subswot_ow_wnl.geojson')
            print('SUCCESS!')


    return riverSP_gdf