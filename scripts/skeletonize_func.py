# FUNCTIONS USED IN SKELETONIZATION PROCESS OF BRAIDED RIVERS 
# (pixel and raster operation functions for skeletonization)

# Import libraries
import os
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, LineString, MultiLineString, Polygon
import numpy as np
import matplotlib.pyplot as plt
import glob
from sklearn.metrics.pairwise import cosine_similarity

from tqdm import tqdm

# Specifically for skeleton pruning
from skimage.morphology import skeletonize, binary_dilation, binary_erosion, label, remove_small_holes
from skimage.filters import gaussian
from scipy.signal import convolve2d
from scipy.spatial import cKDTree
import sys
sys.path.append("/Volumes/OneTouch/work/scripts")
import procBraided as pc

# Importing centerlines and converting between image coordinates and lat lon
import rasterio

# For merging and selecting centerline channels
from sklearn.metrics.pairwise import cosine_similarity
from shapely import centroid
from shapely.ops import linemerge, nearest_points

import tools

def export_progress(array, ref_file, export_file):
    with rasterio.open(ref_file) as src:

        profile = src.profile  # Store metadata
        data = src.read()  # Read all bands

    # overwrite the data
    data[0] = array

    # Export the updated data
    with rasterio.open(export_file, "w", **profile) as dst:
        dst.write(data)



# ---------------------------- (PRIMARY FUNCTIONS) ----------------------------
def get_watermask(path):

    # Extract water mask data from provided geotiff using rasterio
    with rasterio.open(path) as water_mask:
        water_mask_data = water_mask.read(1)  # Read the first band (binary mask)

    return water_mask_data


def smooth_mask(mask_image, pixcdate, figdir, dilate_amount, gauss_amount, gauss_thresh=0.6, max_hole_size=300, season='LF', second_dilation=True, savePlot=False):
    """Function for extracting the initial skeleton from binary raster image
    # INPUTS:
    #       - mask_image: binary raster image of water (1s) and not water (0s)
    #       - pixcdate: string of date 'yyyymmdd' of the pixel cloud acquisition
    #       - figdir: figure directory
    #       - dilate_amount: parameter controlling pixel dilation (meters)
    #       - gauss_amount: parameter controlling gaussian filter kernal size (meters)
    # OUTPUTS:
    #       - skeleton: binary image of returned skeleton
    """

    # create progress bar
    bar = tqdm(total=7, desc="Smoothing mask", leave=False)

    # Step 1: buffer water mask to combine small geometries (across dams, bridges, etc)
    bar.update(1)
    if dilate_amount > 0:
        dilated_mask = binary_dilation(mask_image, footprint=tools.createKernel(dilate_amount)) # OG 10
        dilated_mask = binary_erosion(dilated_mask, footprint=tools.createKernel(dilate_amount)) # erode by same amount
    else:
        dilated_mask = mask_image    
    
    # Step 2: Select only the largest connected structure
    bar.update(1)
    largest_component_mask = tools.find_largest_component(dilated_mask)
    
    # Step 3. Perform gaussian smoothing
    bar.update(1)
    smoothed_mask = gaussian(largest_component_mask.astype(float), sigma=float(gauss_amount)) # OG 20
    smoothed_mask = smoothed_mask > gauss_thresh
    
    # Step 4. Fill small holes in remaining mask
    bar.update(1)
    smoothed_mask = remove_small_holes(smoothed_mask, max_hole_size)

    # Optional Step: Perform secondary dilation
    bar.update(1)
    if second_dilation:
        smoothed_mask = binary_dilation(smoothed_mask, footprint=tools.createKernel(dilate_amount)) # OG 10
        smoothed_mask = binary_erosion(smoothed_mask, footprint=tools.createKernel(dilate_amount)) # erode by same amount
    
    # Step 5: Select again the largest connected structure
    bar.update(1)
    largest_component_mask_smoothed = tools.find_largest_component(smoothed_mask)

    bar.update(1)
    return largest_component_mask_smoothed




def get_skeleton(smoothed_mask, mask_image, pixcdate, figdir,distance_threshold=30, prunethresh=600):

    """
    Skeleton pruning algorithm for multi-channel and braided river centerlines
    
    Parameters:
    skeleton : output skeleton from the get_skeleton function, derived from S2 watermask rasters
    pixcdate : string of format 'yyyymmdd'
    fidir :    Figure directory
    prunethresh : threshold pixel length for pruning branches

    Returns:
    pruned_skeleton : The final pruned skeleton as a np.array() (in image coordinates)

    Method Applied:
    1. Fill small holes based on input hole size
    2. Compute end and joint points (for plotting)
    3. Perform initial branch pruning based on threshold branch length
    4. Iteratively reskeletonize and prune branches until length of endpoints is unchanging
    5. Plot results

    """
    # Create progress bar
    bar = tqdm(total=4, desc="Pruning", leave=False)

    # Step 1. Perform skeletonization
    bar.update(1)
    skeleton = skeletonize(smoothed_mask)

    # Step 2. Connect close endpoints
    bar.update(1)
    connected_skeleton = tools.connect_close_endpoints(skeleton, distance_threshold) # pixels (1 pixel = 10 m)

    # Step 3. Iterate through pruning algortihm until no further changes
    bar.update(1)
    pruned_skeleton = tools.prune_branches(connected_skeleton, prunethresh)
    pruned_skeleton = skeletonize(pruned_skeleton)
    endpoints_len_new = float('inf')

    while True:
        endpoints = tools.find_endpoints(pruned_skeleton)
        y, x = np.where(endpoints == 1) 
        endpoints_len = len(x)

        if endpoints_len_new == endpoints_len:
            final_skeleton = skeletonize(pruned_skeleton) # get rid of areas with 2 pixels width
            break
        
        pruned_skeleton = tools.prune_branches(pruned_skeleton, prunethresh) # APPLY pruning algorithm

        # Compute length of endpoints
        y, x = np.where(endpoints == 1) 
        endpoints_len_new = len(x)

    bar.update(1)
    return final_skeleton


def get_labeled_skeleton(final_skeleton, pixcdate, figdir, savePlot=False):
    
    # Now, set ID values to each channel....
    # endpoints = tools.find_endpoints(final_skeleton)
    joints, joints_removed_mask = tools.find_joints(final_skeleton)
    # id = np.where(endpoints)[0] == np.min(np.where(endpoints)[0])
    # start_pixel = (np.where(endpoints)[0][id][0], np.where(endpoints)[1][id][0]) # (y,x), choose endpoint with the smallest y value (most north endpoint)

    joint_coords = np.where(joints == 1)  # Get skeleton coordinates
    joint_coords = np.dstack((joint_coords[0],joint_coords[1]))
    joint_coords = joint_coords.reshape(-1, 2)

    # Assign labels
    labeled_skeleton = tools.assign_unique_ids_to_branches(final_skeleton, joint_coords, None)

    return labeled_skeleton



def find_nearest_branch(subskel,labeled_skeleton):

    # subsket is a skeleton with less only 1 pixel
    # labeled_skeleton is the original skeleton with branch IDs
    # returns the nearest branch ID to the subskel
    # print('Finding nearest branch to branch with less than 1 pixel...')
    y,x = np.where(subskel > 0)
    y = y[0]
    x = x[0]

    # nearest positive value in labeled_skeleton that is not at the same location as itself, extract branch ID
    # compute distance from y,x to all other postive values in labeled_skeleton
    y2,x2 = np.where(labeled_skeleton > 0)
    distances = np.sqrt((y2-y)**2 + (x2-x)**2)
    distances = np.where(distances == 0, np.inf, distances)
    nearest_branch = np.argmin(distances)

    # print('nearest branch: ',labeled_skeleton[y2[nearest_branch]][x2[nearest_branch]])
    return labeled_skeleton[y2[nearest_branch]][x2[nearest_branch]]



def extract_cl_from_skeleton(labeled_skeleton,water_mask_tiff):

    # print('...........Running centerline extraction from skeleton algorithm...........')

    line_data = []
    skipped_branches = []
    for branch in tqdm(np.unique(labeled_skeleton), desc='Vectorizing product',  leave=False):
        # print('Processing branch No. '+str(branch)+'...')
        if branch == 0:
            continue

        # extract subskeleton
        subskel = np.zeros_like(labeled_skeleton)
        subskel[labeled_skeleton == branch] =  1
        # if subskeleton has less than 1 pixels, store the branch ID and append that pixel to nearest branch outside of the branch loop
        if len(np.where(subskel > 0)[0]) == 1:
            # print('Branch No. '+str(branch)+', has only than 1 pixel! Skipping and merging later...')
            skipped_branches.append(branch)
            continue

        # order raster from endpoint
        pixel_coords_x, pixel_coords_y = tools.order_raster_from_endpoint(subskel,branch)
        if pixel_coords_x is None:
            continue

        # Convert ordered pixel coordinates to lat lon
        line = tools.pixel_coordinates2latlonline(water_mask_tiff,pixel_coords_y,pixel_coords_x,branch,reverse=True)

        # append gdf line to overall gdf with branch id
        line_data.append({'branch_id': branch, 'geometry': line})

    # go through skipped branches, and append to nearest branch
    for branch in tqdm(skipped_branches, desc='Vectorizing skipped branches',  leave=False):
        # print('Merging pixels in branch No. '+str(branch)+' to nearest branch...')
        # extract subskeleton
        subskel = np.zeros_like(labeled_skeleton)
        subskel[labeled_skeleton == branch] =  1

        # locate nearest pixel from another branch
        nearest_branch = find_nearest_branch(subskel,labeled_skeleton)

        # merge pixels
        labeled_skeleton[labeled_skeleton == branch] = nearest_branch

        # rerun order raster from endpoint and convert to lat lon for the nearest branch with new appended pixel
        subskel = np.zeros_like(labeled_skeleton)
        subskel[labeled_skeleton == nearest_branch] =  1
        pixel_coords_x,pixel_coords_y = tools.order_raster_from_endpoint(subskel,nearest_branch)

        # Convert ordered pixel coordinates to lat lon
        line = tools.pixel_coordinates2latlonline(water_mask_tiff,pixel_coords_y,pixel_coords_x,nearest_branch,reverse=True)

        # REMOVE old line with branch_id = nearest_branch from line_data
        line_data = [line for line in line_data if line['branch_id'] != nearest_branch]

        # append gdf line to overall gdf with branch id
        line_data.append({'branch_id': nearest_branch, 'geometry': line})

    cl_gdf_labeled = gpd.GeoDataFrame(line_data, crs="EPSG:4326")  # Assuming WGS84
    return cl_gdf_labeled



def extract_cl_from_skeleton_old(labeled_skeleton,water_mask_tiff):

    # print('...........Running centerline extraction from skeleton algorithm...........')

    line_data = []
    for branch in np.unique(labeled_skeleton):
        if branch == 0:
            continue

        # extract subskeleton
        subskel = np.zeros_like(labeled_skeleton)
        subskel[labeled_skeleton == branch] =  1

        # order raster from endpoint
        pixel_coords_x,pixel_coords_y = tools.order_raster_from_endpoint(subskel,branch)

        # Convert ordered pixel coordinates to lat lon
        line = tools.pixel_coordinates2latlonline(water_mask_tiff,pixel_coords_y,pixel_coords_x,branch,reverse=True)

        # append gdf line to overall gdf with branch id
        line_data.append({'branch_id': branch, 'geometry': line})

    cl_gdf_labeled = gpd.GeoDataFrame(line_data, crs="EPSG:4326")  # Assuming WGS84
    #cl_gdf_labeled.plot()
    return cl_gdf_labeled

def calculate_vector(branch_geom, endpoint_index):
    # Calculate the direction vector for a LineString at the given endpoint.
    coords = list(branch_geom.coords)
    if endpoint_index == 0:  # Start of the branch
        vec = np.array(coords[1]) - np.array(coords[0])
    else:  # End of the branch
        vec = np.array(coords[-1]) - np.array(coords[-2])
    return vec / np.linalg.norm(vec)  # Normalize


def get_endpoints_and_vectors(linestring):
    """Return the endpoints and their respective vectors."""
    return [
        (tuple(linestring.coords[0]), calculate_vector(linestring, 0)),  # Start point
        (tuple(linestring.coords[-1]), calculate_vector(linestring, -1)),  # End point
    ]


def merge_short_centerlines(gen_cl, hemi, connection_threshold=50, similarity_threshold=0.7):
    # Created with help from CHAT GPT
    # Prompt:
    # 1. function input is the geopandas gdf with each entry having a unique branch id and linestring geometry
    # 2. convert gdf to utm coordinate system
    # 3. Identify all endpoints for each branch
    # 4. Identify all branches with branch length < 2000 m
    # 5. Go through a merge branching routine, where starting from the first 'short' branch, we identify all branch endpoints within 
    # a threshold distance of each endpoint of the current branch (so that we only search for connecting branches) and we merge the small 
    # branch with the connecting branch that has also the greatest directional similarity (by checking the cosine similarity of their vectors). 
    # - However, if there exists no directional similarity greater than at least 0.7, we leave the current branch as is and do not merge 
    # it with any other branch. Then we move on, even though we left the last branch is left with a length < 2000 m. 
    # - If on the other hand we do merge the short branch with a connection that also has directional similarity > 0.7, we need to next 
    # recompute the list of short branches, now considering that we have made the merge. In this way, we forget the original geodataframe 
    # (the old one where nothing has been merged) and continue forward with the new updated one. This needs to be continued until the 
    # geodataframe is unchanging (no more short branches, that also have a directional similarity > 0.7 with another branch, are identified.)

    # Merges branches based on directional similarity of endpoints (similarity_threshold) and proximity from endpoint to endoint (connection_threshold)

    # 1. Convert to UTM CRS
    centroid = gen_cl.geometry.unary_union.centroid
    utm_zone = int((centroid.x + 180) // 6) + 1
    if hemi == 'north':
        utm_crs = f'EPSG:{32600 + utm_zone}'
    elif hemi == 'south':
        utm_crs = f'EPSG:{32700 + utm_zone}'
    else:
        raise ValueError("Invalid hemisphere. Use 'north' or 'south'.")
    
    cl_utm = gen_cl.to_crs(utm_crs)
    cl_utm['length'] = cl_utm.length

    # 2. Iteratively merge branches
    while True:
        short_branches = cl_utm[cl_utm['length'] < 2000]
        if short_branches.empty:
            break  # Stop if no short branches remain

        for branch_id, branch_geom in short_branches[['branch_id', 'geometry']].values:
            branch_endpoints = get_endpoints_and_vectors(branch_geom)
            closest_candidates = []

            # Find connecting branches within the threshold distance
            for candidate_id, candidate_geom in cl_utm[['branch_id', 'geometry']].values:
                if branch_id == candidate_id:
                    continue

                candidate_endpoints = get_endpoints_and_vectors(candidate_geom)

                for branch_endpoint, branch_vector in branch_endpoints:
                    for candidate_endpoint, candidate_vector in candidate_endpoints:
                        distance = Point(branch_endpoint).distance(Point(candidate_endpoint))
                        if distance <= connection_threshold:
                            similarity = cosine_similarity(
                                branch_vector.reshape(1, -1), 
                                candidate_vector.reshape(1, -1)
                            )[0, 0]
                            closest_candidates.append({
                                "branch_id": candidate_id,
                                "distance": distance,
                                "similarity": similarity,
                                "branch_endpoint": branch_endpoint,
                                "candidate_endpoint": candidate_endpoint
                            })

            # Select the best candidate based on similarity
            if closest_candidates:
                best_candidate = max(
                    closest_candidates,
                    key=lambda x: (x['similarity'], -x['distance'])  # Prioritize similarity, then proximity
                )
                if best_candidate['similarity'] > similarity_threshold:
                    # Merge branches
                    candidate_id = best_candidate["branch_id"]
                    candidate_geom = cl_utm.loc[cl_utm['branch_id'] == candidate_id, 'geometry'].iloc[0]

                    # Combine coordinates
                    all_coords = []
                    if best_candidate['branch_endpoint'] == tuple(branch_geom.coords[0]):
                        all_coords.extend(list(branch_geom.coords[::-1]))
                    else:
                        all_coords.extend(list(branch_geom.coords))

                    if best_candidate['candidate_endpoint'] == tuple(candidate_geom.coords[0]):
                        all_coords.extend(list(candidate_geom.coords))
                    else:
                        all_coords.extend(list(candidate_geom.coords[::-1]))

                    # Update GeoDataFrame
                    merged_geom = LineString(all_coords)
                    cl_utm = cl_utm[cl_utm['branch_id'] != branch_id]
                    cl_utm = cl_utm[cl_utm['branch_id'] != candidate_id]

                    new_row = pd.DataFrame({
                        "branch_id": [candidate_id],  # Retain candidate ID (the larger branch)
                        "geometry": [merged_geom],
                        "length": [merged_geom.length]
                    })
                    cl_utm = pd.concat([cl_utm, gpd.GeoDataFrame(new_row, crs=cl_utm.crs)])
                    break  # Restart the loop after updating

        else:
            # If no merges occurred, exit the loop
            break

    cl_new  = cl_utm.sort_values(by=['branch_id'])
    cl_new = cl_new.to_crs(gen_cl.crs)
    cl_new = cl_new.reset_index(drop=True)
    cl_new['branch_id_old'] = cl_new['branch_id']
    cl_new['branch_id'] = cl_new.index + 1


    # Return the updated GeoDataFrame
    return cl_new



def trim_cl_to_river_bounds(cl, river_bounds):
    """Depending on the watermask and starting point of the river system, there can be circular loops generated at the start of the system.
    By trimming the centerline to a start line (river cross section) we can ensure that the starting reach(es) is(are) clean.

    Parameters
    ----------
    cl : _type_
        _description_
    starting_line : _type_
        _description_
    """

    # cut the geodataframe by the centerline
    cl_new = cl.copy()
    # cl_new['geometry'] = cl.within(river_bounds)
    cl_new = gpd.clip(cl_new, river_bounds)

    return cl_new


from shapely.ops import polygonize, unary_union, linemerge
from shapely import MultiPolygon, MultiPoint, concave_hull
import time

def join_cl_at_joints(cl, starting_line, main_centerline, search_dist=20):

    """Reaches are seperated due to the processing of pixels and selecting their center coordinates. we join them here and also keep track of which reach each reach flows into"""

    # fist convert centerline to local coordinates
    local_crs = cl.estimate_utm_crs()
    cl = cl.to_crs(local_crs)
    starting_line = starting_line.to_crs(local_crs)
    main_centerline = main_centerline.to_crs(local_crs)
    # main_centerline = linemerge(main_centerline.unary_union)
    main_centerline = main_centerline.geometry[0]


    # see if we need to reverse coordinates in main centerline
    main_centerline_endpoints = main_centerline.coords[0], main_centerline.coords[-1]
    dist = [Point(coord).distance(starting_line.unary_union) for coord in main_centerline_endpoints]
    if dist[0] > dist[-1]: # reverse
        reversed_coords = list(main_centerline.coords)
        reversed_coords.reverse()
        main_centerline = LineString(reversed_coords)

    
    # Cycle through centerlines, and document endpoints
    branch_inlet = list()
    branch_outlet = list()
    branch_proj_dist = list()
    for i in cl.index:

        branch_id = cl.loc[i, 'branch_id']
        line = cl.loc[i, 'geometry']

        # extract coords, and order from first to last according to main center line
        dist   = [main_centerline.project( Point(coord), normalized=True ) for coord in [line.coords[0], line.coords[-1]]]

        # see if we need to reverse the line string
        if dist[0] > dist[-1]: # reverse
            reversed_coords = list(line.coords)
            reversed_coords.reverse()
            line = LineString(reversed_coords)
            cl.loc[i, 'geometry'] = line

        # save endpoints and assign the reach start dist and the reach dist
        inlet = Point(line.coords[0])
        outlet = Point(line.coords[-1]) # order the points
        proj_dist = np.min(dist)

        # add values to list
        branch_inlet.append(inlet)
        branch_outlet.append(outlet)
        branch_proj_dist.append(proj_dist)

    # add to dataframe, create new index and sort by distance
    inlet_coords = dict(zip(cl.branch_id.values, branch_inlet))
    outlet_coords = dict(zip(cl.branch_id.values, branch_outlet))
    cl['branch_proj_dist'] = branch_proj_dist
    cl = cl.sort_values(by='branch_proj_dist').reset_index(drop=True)

    # now cycle through dataframe and assign downstream ids, only assigning to upstream points
    cl['downstream_ids'] = None
    for i in cl.index:

        branch_id = cl.loc[i, 'branch_id']

        # Check what inlet coords are within the search distance if the branches outlet
        downstream_branches = list()
        for key in inlet_coords.keys():
            if key != branch_id:
                if outlet_coords[branch_id].distance(inlet_coords[key]) < search_dist:
                    downstream_branches.append(key)

        # Check also if there are neighboring outlet branches that should be joined
        shared_outlet_branches = list()
        for key in inlet_coords.keys():
            if key != branch_id:
                if outlet_coords[branch_id].distance(outlet_coords[key]) < search_dist:
                    shared_outlet_branches.append(key)

        # add the list of downstram brach ids to the dataframe
        cl.loc[i, 'downstream_ids'] = ', '.join([str(key) for key in downstream_branches])



        # Find the new common jon for the outlet and downstream inlets
        shared_point = MultiPoint([inlet_coords[key] for key in downstream_branches]+ [outlet_coords[key] for key in shared_outlet_branches] + [outlet_coords[branch_id]]).centroid
        outlet_coords[branch_id] = shared_point
        for key in downstream_branches:
            inlet_coords[key] = shared_point
        for key in shared_outlet_branches:
            outlet_coords[key] = shared_point

        
        
    # Finally, update endpoints to the revised editions
    # loop through each branch, and update the endpoints with whats in the updated dictionaries
    for i in cl.index:

        # extract points in linestring for branch
        branch_id = cl.loc[i, 'branch_id']
        coords = [Point(coord) for coord in cl.loc[i, 'geometry'].coords]
        coords[0] = inlet_coords[branch_id]
        coords[-1] = outlet_coords[branch_id]

        # update in place
        cl.loc[i, 'geometry'] = LineString(coords)


    return cl

    # Find the starting reaches 