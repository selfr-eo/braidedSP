#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Functions designed for processing and working with SWOT products over rivers and lakes

from datetime import datetime as dt
import geopandas as gpd
import pandas as pd
from shapely.geometry import LineString, Point
#from sklearn.cluster import DBSCAN
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import glob as glob
#import netCDF4 as nc
# For cl ordering algorithm
import networkx as nx
from shapely import centroid

from scipy.spatial import cKDTree

# For trimming to raster mask
import rasterio
from rasterio.transform import rowcol




# ---------------------------- centerline ordering functions ----------------------------


def haversine(lat1, lon1, lat2, lon2):
    import math
    """
    Returns:
    distance : Distance between the two points in meters
    """
    
    # Radius of Earth in meters
    R = 6371000.0
    
    # Convert latitude and longitude from degrees to radians
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)
    
    # Difference in coordinates
    dlat = lat2_rad - lat1_rad
    dlon = lon2_rad - lon1_rad
    
    # Haversine formula
    a = math.sin(dlat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    
    # Distance in meters
    distance = R * c
    
    return distance

def sort_SWORD_cl(cl,showPlots=True):

    # 1. EXTRACT REACH ENDPOINTS

    line_endpoints = []
    for geom in cl.geometry:
        # ASSUMING coordinates are listed in order from start to end for each line...
        coords = list(geom.coords)
        start_point = coords[0]
        end_point = coords[-1]
        line_endpoints.append([start_point, end_point])  # Store (start, end)

    # 2. COMPUTE DISTANCE FROM EACH ENDPOINT TO ALL OTHER ENDPOINTS

    distance_list = []
    pixelID = [] # used for tracking endpoints

    for geomID in range(len(cl.geometry)):

        # extract start and endpoint coordinates for current reach (geom)
        slon, slat = line_endpoints[geomID][0]
        elon, elat = line_endpoints[geomID][1]

        # set values for closest distance to start, end, and closest index
        ds = float('inf')
        de = float('inf')
        dmin = float('inf')
        dmax = float(0)

        connectingColID = 1     # Assume the current reach connects via the end point unless said otherwise
        #closestColID = 0        # Take closest pixel is a start pixel unless said otherwise
        #closestReachID = -1     # Set index of closest reach
        closestColID_s = closestColID_e = closestReachID_s = closestReachID_e = -1

        for i in range(len(cl.geometry)):

            if i == geomID:
                continue # skip distance calculation for itself

            slon2, slat2 = line_endpoints[i][0] # start point of next geom
            elon2, elat2 = line_endpoints[i][1] # endpoint of next geom

            dss_se_ee_es = [haversine(slat,slon,slat2,slon2), haversine(slat,slon,elat2,elon2),haversine(elat,elon,elat2,elon2),haversine(elat,elon,slat2,slon2)]

            distance = min(dss_se_ee_es)
            ds_min = min(dss_se_ee_es[0:2])
            de_min = min(dss_se_ee_es[2:4])

            # to identify endpoints: start and end will be connected to same point when a reach has an endpoint
            if ds_min < ds:
                ds = ds_min
                closestReachID_s = i
                closestColID_s = 0 if dss_se_ee_es.index(ds_min) == 0 else 1
            if de_min < de:
                de = de_min
                closestReachID_e = i
                closestColID_e = 0 if dss_se_ee_es.index(de_min) == 3 else 1

            candidate_ep = 1 if (closestColID_s == closestColID_e) and (closestReachID_s == closestReachID_e) and (closestReachID_s != -1) and (closestReachID_e != -1) else 0 # Mark as candidate endpoint

            if distance < dmin:
                # Find absolute closest connection

                dmin = distance
                min_index = dss_se_ee_es.index(dmin)

                if min_index < 2:
                    # Min distance of current reach is connected via the current start pixel
                    connectingColID = 0
                else:
                    # Min distance of current reach is connected via the current end pixel
                    connectingColID = 1

                closestColID = 0 if min_index == 0 or min_index == 3 else 1
                closestReachID = i

        # Append data to distance list
        distance_list.append([geomID,ds,de,closestReachID_s,closestReachID_e,closestColID_s,closestColID_e,connectingColID,candidate_ep])

    # Find true endpoints of entire CL
    true_endpoint_reaches = [entry for entry in distance_list if entry[8] == 1]

    if len(true_endpoint_reaches) != 2:
        raise ValueError("Centerline error: More or less than 2 endpoints found")

    if showPlots == True:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(cl.get_coordinates().x,cl.get_coordinates().y,c='blue',s=0.5)
        for i in range(len(true_endpoint_reaches)):

            if true_endpoint_reaches[i][7] == 0: # if connectingColID == 0 (point with smallest distance is start), endpoint = 1 (end)
                ep1 = line_endpoints[true_endpoint_reaches[i][0]][1]
            if true_endpoint_reaches[i][7] == 1: # if connectingColID == 1 (point with smallest distance is the end), endpoint = 0 (start)
                ep1 = line_endpoints[true_endpoint_reaches[i][0]][0]
            if i == 0:
                ax.plot(ep1[0],ep1[1],'og')
            if i == 1:
                ax.plot(ep1[0],ep1[1],'or')
                ax.set_title('Found endpoints')
        plt.show()

    # SELECT FIRST ENDPOINT 
    if true_endpoint_reaches[0][7] == 0: # if connectingColID == 0 (point with smallest distance is start), endpoint = 1 (end)
            true_endpoint = line_endpoints[true_endpoint_reaches[0][0]][1]
    if true_endpoint_reaches[0][7] == 1: # if connectingColID == 1 (point with smallest distance is the end), endpoint = 0 (start)
            true_endpoint = line_endpoints[true_endpoint_reaches[0][0]][0]

    sorted_coords, reach_order_id = sort_cl_by_closestReach(true_endpoint_reaches[0][0],true_endpoint_reaches[1][0],distance_list,cl)

    lon_coordinates = [coord[0] for sublist in sorted_coords for coord in sublist]
    lat_coordinates = [coord[1] for sublist in sorted_coords for coord in sublist]
    if showPlots == True:
        plt.figure()
        plt.plot(lon_coordinates,lat_coordinates)
        plt.title('Sorted centerline coordinates')


    cl_gdf = gpd.GeoDataFrame()
    cl_gdf_points = gpd.GeoDataFrame()

    i = 0

    # Save as gdf with lines, also with points
    for sublist in sorted_coords:
        lon_coordinates = [coord[0] for coord in sublist]
        lat_coordinates = [coord[1] for coord in sublist]

        #df_temp = pd.DataFrame({'lon':lon_coordinates,'lat':lat_coordinates})
        coordinates = list(zip(lon_coordinates, lat_coordinates))
        line = LineString(coordinates)
        gdf_temp = gpd.GeoDataFrame(geometry=[line],crs='EPSG:4326')

        # Saving the centroid only
        points_temp = gpd.GeoDataFrame(geometry=[centroid(line)],crs='EPSG:4326')
        points_temp['reachID'] = cl.reach_id[reach_order_id[i]]
        points_temp['width'] = cl.width[reach_order_id[i]]


        gdf_temp['reachID'] = cl.reach_id[reach_order_id[i]]
        gdf_temp['width'] = cl.width[reach_order_id[i]]

        i+=1

        cl_gdf = pd.concat([cl_gdf,gdf_temp])
        cl_gdf_points = pd.concat([cl_gdf_points,points_temp])

    # save as a line
    coordinates = list(zip(cl_gdf.get_coordinates().x, cl_gdf.get_coordinates().y))
    line = LineString(coordinates)

    reverse = False

    if reverse == True:
        coords_df = cl_gdf.get_coordinates()
        coordinates = list(zip(coords_df['x'].iloc[::-1], coords_df['y'].iloc[::-1]))
        line = LineString(coordinates)
        cl_gdf = gpd.GeoDataFrame(geometry=[line],crs='EPSG:4326')

    return cl_gdf, line, cl_gdf_points




# ---------------------------- other data processing functions ----------------------------

def getNearestMaskDate(pixcdate, tileID, LFmonths, HFmonths):
    # get nearest river mask date for a given pixel cloud data date
    # inputs:
    # pixcdate - date of pixel cloud data acquisition in the format 'yyyymmdd'
    # tileID - granule id in format 'PASS_TILE' ex 258_112L
    # LFmonths and HFmonths - list of low flow and high flow months, respectively

    maskdate = pixcdate[4:6]+pixcdate[2:4]
    water_mask_tiff = glob.glob("/Volumes/OneTouch/work/water_masks/brahmaputra/"+tileID+"/S2*"+str(maskdate)+".tif")


     # Need to find next closest mask in the same season...
    if not water_mask_tiff:
        currentmonth = pixcdate[4:6]
        if int(currentmonth) in LFmonths:
            season_months = LFmonths
        else:
            season_months = HFmonths


        dist = np.abs(np.array(season_months) - int(currentmonth))
        for j in np.unique(dist):
            if dist[dist == j][0] == 0: # skip if comparing with self
                continue
            if len(dist[dist == j]) > 1:
                # check first closest month
                monthid_indist = np.where(dist == j)[0][0]
                nextmonth = season_months[monthid_indist]
                if nextmonth < 10:
                    nextmonth = '0'+str(nextmonth)
                else:
                    nextmonth = str(nextmonth)

                maskdate = nextmonth+pixcdate[2:4]
                water_mask_tiff = glob.glob("/Volumes/OneTouch/work/water_masks/brahmaputra/"+tileID+"/S2*"+str(maskdate)+".tif")
                if water_mask_tiff:
                    break

                # check second closest month
                monthid_indist = np.where(dist == j)[0][1]
                nextmonth = season_months[monthid_indist]
                if nextmonth < 10:
                    nextmonth = '0'+str(nextmonth)
                else:
                    nextmonth = str(nextmonth)

                maskdate = nextmonth+pixcdate[2:4]
                water_mask_tiff = glob.glob("/Volumes/OneTouch/work/water_masks/brahmaputra/"+tileID+"/S2*"+str(maskdate)+".tif")
                if water_mask_tiff:
                    break

            else:
                monthid_indist = np.where(dist == j)[0][0]
                nextmonth = season_months[monthid_indist]
                if nextmonth < 10:
                    nextmonth = '0'+str(nextmonth)
                else:
                    nextmonth = str(nextmonth)

                maskdate = nextmonth+pixcdate[2:4]
                water_mask_tiff = glob.glob("/Volumes/OneTouch/work/water_masks/brahmaputra/"+tileID+"/S2*"+str(maskdate)+".tif")
                if water_mask_tiff:
                    break

    print('Closest mask date:',str(maskdate))
    return maskdate, water_mask_tiff[0]


def projectToCenterline(gdf_cl, gdf_dat, hemi):

    # Determine utm zone
    centroid = gdf_cl.geometry.unary_union.centroid
    utm_zone = int((centroid.x + 180) // 6) + 1  # Calculate UTM zone
    if hemi == 'north':
        utm_crs = f'EPSG:{32600 + utm_zone}'
    if hemi == 'south':
        utm_crs = f'EPSG:{32700 + utm_zone}'

    # Reproject to UTM
    gdf_cl = gdf_cl.to_crs(utm_crs)
    gdf_dat_utm = gdf_dat.to_crs(utm_crs)


    # Extract points
    datutm = gdf_dat_utm.get_coordinates()
    clutm = gdf_cl.get_coordinates()

    #centerline = gdf_cl.geometry
    line = LineString(clutm)
    datutm['point_geom'] = datutm.apply(lambda row: Point(row['x'], row['y']), axis=1)
    dist = datutm['point_geom'].apply(line.project)
    return dist 


def readPIXC(filename,gdf_buffered):

    # If no centerline buffer to trim to, set gdf_buffered = []
    
    nc = xr.open_mfdataset(filename, group = 'pixel_cloud', engine='h5netcdf')
    
    
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


    gdf = gpd.GeoDataFrame(pd.DataFrame({"height":height,"heightEGM":heightEGM,"geoid":geoid,"lat":lat_flat,"lon":lon_flat,"class":class_flat[class_flat == 4],"water_frac":water_frac,"phase_noise_std":phase_noise_std,"dheight_dphase":dheight_dphase,"sig0":sig0}),geometry=gpd.points_from_xy(lon_flat,lat_flat))
    gdf.set_crs(epsg=4326, inplace=True)


    ## Trim data to buffer around imported centerline
    if not gdf_buffered:
        clipped_gdf = gdf
        print("No centerline buffer to trim to.")
    else:
    #if gdf_buffered.empty == False:
        print('Trimming data to centerline buffer..')
        clipped_gdf = gpd.clip(gdf, gdf_buffered)

    return clipped_gdf


def trim2mask_general(pixc_gdf, mask_nparray, mask_tiff_filename):
    # This function created with help from ChatGPT
    # - Function to trim pixel cloud data to a binary mask (numpy array of 1s and 0s)
    # - Water mask generated from S2 optical imagery and imported as a geotiff
    # Open the GeoTIFF file
    with rasterio.open(mask_tiff_filename) as water_mask:
        mask_transform = water_mask.transform  # Get the affine transform
        mask_crs = water_mask.crs  # Get the CRS (coordinate reference system)

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
        if 0 <= row < mask_nparray.shape[0] and 0 <= col < mask_nparray.shape[1]:
            sampled_values.append(mask_nparray[row, col])
        else:
            sampled_values.append(0) 

    pixc_gdf['water_masked'] = np.array(sampled_values) == 1

    # Filter the GeoDataFrame to keep only points within the mask
    trimmed_gdf = pixc_gdf[pixc_gdf['water_masked']]

    # Plot masked results
    pixc_gdf.plot(
        column='water_masked',
        markersize=0.1,
        legend=True,
        categorical=True,
        cmap='tab20'
    )
    
    return trimmed_gdf, pixc_gdf