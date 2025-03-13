#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# FUNCTIONS USED IN SKELETONIZATION PROCESS OF BRAIDED RIVERS 
# (pixel and raster operation functions for skeletonization)

# Import libraries
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, LineString, MultiLineString
import numpy as np
import matplotlib.pyplot as plt
import glob
from sklearn.metrics.pairwise import cosine_similarity

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



# ---------------------------- (HELPER FUNCTIONS) ----------------------------
# (Functions in order of appearance)


def get_watermask(water_mask_tiff):
    # Extract water mask data from provided geotiff using rasterio

    #import rasterio

    with rasterio.open(water_mask_tiff) as water_mask:
        water_mask_data = water_mask.read(1)  # Read the first band (binary mask)
        # mask_transform = water_mask.transform  # Get the affine transform
        # mask_crs = water_mask.crs  # Get the CRS (coordinate reference system)
        # season = 'LF'

    return water_mask_data


def createKernel(radius):
    # create circular kernel for gaussian smoothing
    kernel = np.zeros((2*radius+1, 2*radius+1))
    y,x = np.ogrid[-radius:radius+1, -radius:radius+1]
    mask = x**2 + y**2 <= radius**2
    kernel[mask] = 1
    return kernel


def find_joints(skeleton):
    # (do not count the center pixel)
    kernel = np.array([[1, 1, 1],
                       [1, 0, 1],
                       [1, 1, 1]])

    neighbor_count = convolve2d(skeleton.astype(int), kernel, mode='same')
    jointpoints = (neighbor_count > 2) & (skeleton.astype(int) == 1)
    joints_removed_mask = (skeleton.astype(int) == 1) & (neighbor_count <= 2)

    return jointpoints, joints_removed_mask

def find_endpoints(skeleton):
    # (do not count the center pixel)
    kernel = np.array([[1, 1, 1],
                       [1, 0, 1],
                       [1, 1, 1]])

    neighbor_count = convolve2d(skeleton.astype(int), kernel, mode='same')
    endpoints = (neighbor_count == 1) & (skeleton.astype(int) == 1)

    return endpoints


def find_indirect_neighbors(skeleton):
    # (do not count the center pixel)
    kernel = np.array([[1, 0, 1],
                       [0, 0, 0],
                       [1, 0, 1]])

    neighbor_count = convolve2d(skeleton.astype(int), kernel, mode='same')
    endpoint_candidates = (neighbor_count == 1) & (skeleton.astype(int) == 1)

    return endpoint_candidates


def trace_branch_from_endpoint(endpoint,joint_coords,pruned_skeleton):
    # Function to trace branch starting from an endpoint

    branch_pixels = [endpoint]
    visited = set(branch_pixels)
    current_pixel = endpoint

    while True:
        neighbors = [(current_pixel[0] + i, current_pixel[1] + j) 
                    for i in [-1, 0, 1] for j in [-1, 0, 1] 
                    if not (i == 0 and j == 0)]
        
        # Remove neighbors outside of image
        img_yx_shape = np.shape(pruned_skeleton)
        update_neighbors = [p for p in neighbors if p[0] < img_yx_shape[0] and p[1] < img_yx_shape[1]]

        next_pixels = [p for p in update_neighbors if pruned_skeleton[p] == 1 and p not in visited]
        #next_pixels = [p for p in neighbors if pruned_skeleton[p] == 1 and p not in visited]
            
        if not next_pixels:
            break
            
        # Check if we have reached a joint
        joint_coords_set = {tuple(coord) for coord in joint_coords}
        if any(pixel in joint_coords_set for pixel in next_pixels):
            #print('Pixel is a joint!')
            break
        
        current_pixel = next_pixels[0]
        branch_pixels.append(current_pixel)
        visited.add(current_pixel)

    return branch_pixels


def add_line(skeleton, start, end):
    # Function to add a straight line of skeleton pixels between two points
    y0, x0 = start
    y1, x1 = end
    length = max(abs(y1 - y0), abs(x1 - x0)) + 1
    y_values = np.linspace(y0, y1, length).round().astype(int)
    x_values = np.linspace(x0, x1, length).round().astype(int)
    skeleton[y_values, x_values] = 1
    return skeleton


def connect_close_endpoints(skeleton, distance_threshold=10):
    # Connect close endpoints
    endpoints = find_endpoints(skeleton)
    endpoint_coords = np.column_stack(np.where(endpoints))
    
    if len(endpoint_coords) < 2:
        return skeleton  # No connections possible with fewer than 2 endpoints

    # Use KDTree to find endpoint pairs within the distance threshold
    tree = cKDTree(endpoint_coords)
    pairs = tree.query_pairs(distance_threshold)
    
    for (idx1, idx2) in pairs:
        point1 = tuple(endpoint_coords[idx1])
        point2 = tuple(endpoint_coords[idx2])
        skeleton = add_line(skeleton, point1, point2)

    return skeleton


def fill_small_holes(skeleton, max_hole_size=3):
    # Close small gaps in the skeleton to fill tiny loops
    filled_skeleton = remove_small_holes(skeleton, area_threshold=max_hole_size**2)

    # Fill small holes specifically up to max_hole_size (in pixels)    
    return filled_skeleton.astype(int)


def pixel_coordinates2latlonline(geotiff_path, pixel_y, pixel_x,branch_id,reverse=False):
    # - Function to convert image coordinates to latitude and longitude
    # - pixel_y and pixel_x are already sorted as ascending from their endpoint (first entry)
    # - Returns line geometry of the pixel coordinates converted to latitude and longitude

    pixel_coords = np.column_stack((pixel_y, pixel_x))

    with rasterio.open(geotiff_path) as dataset:
        latlon_coords = np.array([dataset.xy(row, col) for row, col in pixel_coords])

    df = pd.DataFrame({
        'latitude': latlon_coords[:, 1],  # latitude is the second item in (x, y) output from dataset.xy
        'longitude': latlon_coords[:, 0], # longitude is the first item in (x, y)
    })

    points = [tuple(x) for x in df[['longitude', 'latitude']].to_numpy()]

    if len(points) > 1:  # Only create LineString if there are at least two points
        line = LineString(points)
        if reverse == True:
            line = LineString(line.coords[::-1])
    if len(points) <= 1:
        raise ValueError("1 or less points in branch No. {}".format(str(branch_id)))
    
    return line


def latlon_to_pixel_coordinates(geotiff_path, latitudes, longitudes):
    # - Function to convert SWORD lat lon cl to image coordinates
    # - Performs the opposite of function extract_latlon_with_branch_id

    with rasterio.open(geotiff_path) as dataset:
        transform = dataset.transform

    # Convert lat/lon to pixel (row, col) using the inverse of the transform
    pixel_coords = [~transform * (lon, lat) for lat, lon in zip(latitudes, longitudes)]
    
    # Convert the result to a pandas df
    pixel_coords_df = pd.DataFrame(pixel_coords, columns=['x', 'y'])


    
    return pixel_coords_df

# ---------------------------- (INTERMEDIATE FUNCTIONS) ----------------------------

def compute_branchlen(endpoint,joint_coords,pruned_skeleton):
    distances = np.sqrt((joint_coords[:, 0] - endpoint[0])**2 + (joint_coords[:, 1] - endpoint[1])**2)
    # Find the minimum distance
    min_distance = np.min(distances) # in pixels
    return(min_distance)

def prune_branches(skeleton, length_threshold):
    # Function for pruning branches
    endpoints = find_endpoints(skeleton)
    joints, joints_removed_mask = find_joints(skeleton)
    joint_coords = np.where(joints == 1)  # Get skeleton coordinates
    joint_coords = np.dstack((joint_coords[0],joint_coords[1]))
    joint_coords = joint_coords.reshape(-1, 2)
    pruned_skeleton = skeleton.copy()
    miny = 100
    yrange = np.shape(skeleton)[0]
    maxy = yrange - 100

    # Iterate over all endpoints and trace branches
    for i in range(len(np.where(endpoints)[0])):
        endpoint = (np.where(endpoints)[0][i], np.where(endpoints)[1][i]) # y, x
        
        # If endpoint crosses frame exit, keep branch (FOR N-S traveling rivers..need to add option for E-W rivers)
        # If endpoint is in the top or bottom 100 m of the frame, keep!
        # if endpoint[0] == np.min(np.where(endpoints)[0]) or endpoint[0] == np.max(np.where(endpoints)[0]):
        #     continue
        if endpoint[0] < miny or endpoint[0] > maxy:
            continue

        # Only select branches to search over that are less than a threshold length
        branchlen = compute_branchlen(endpoint, joint_coords, pruned_skeleton)
        if branchlen < length_threshold:
            # need to now locate all pixels in this branch and remove them
            branch = trace_branch_from_endpoint(endpoint,joint_coords,pruned_skeleton) #returns branch pixels
            if branch == None:
                print('No pixels in branch!')
                continue
            for pixel in branch:
                pruned_skeleton[pixel] = 0  # Remove branch pixels

    return pruned_skeleton


def trace_branch_and_label(startpoint, initial_joint, joint_coords, pruned_skeleton, labeled_segments, assign_id, visited_segments):
    # This function for tracing !ONE! branch from the starting joint until it reaches a joint or end point
    # The pixel 'startpoint' is one of the neighbors of the initial joint
    # The list 'initial_joint' prevents the branch from going in the direction of the initial joint
    branch_pixels = [startpoint]
    visited = set(branch_pixels)
    initial_joint = set([initial_joint])
    current_pixel = startpoint
    joint_coords_set = {tuple(coord) for coord in joint_coords}

    while True:
        neighbors = [(current_pixel[0] + i, current_pixel[1] + j) 
                    for i in [-1, 0, 1] for j in [-1, 0, 1] 
                    if not (i == 0 and j == 0)]

        # Remove neighbors outside of image
        img_yx_shape = np.shape(pruned_skeleton)
        update_neighbors = [p for p in neighbors if p[0] < img_yx_shape[0] and p[1] < img_yx_shape[1]]

        next_pixels = [p for p in update_neighbors if pruned_skeleton[p] == 1 and p not in visited and p not in initial_joint]
        
        #next_pixels = [p for p in neighbors if pruned_skeleton[p] == 1 and p not in visited and p not in initial_joint]
        
        if not next_pixels:
            break


        if any(pixel in joint_coords_set for pixel in next_pixels):
            break

        current_pixel = next_pixels[0]     # Should just be one next pixel if not a joint
        branch_pixels.append(current_pixel)
        visited.add(current_pixel)

    for pixel in branch_pixels:
        labeled_segments[pixel] = assign_id
        visited_segments.add(pixel)         # Mark as processed

    return labeled_segments, visited_segments


def assign_unique_ids_to_branches(pruned_skeleton, joint_coords, starting_pixel=None): # USE THIS ONE
    ### Wrapper function to assign unique IDs to all branches
    #  -- Starting pixel should be the northern most jointpoint (if provided)
    #  -- Joint_coords are sorted by y index (smallest first)
    #  RETURNS:
    # updated_labeled_segments: a numpy array same size as pruned skeleton with unique branch IDs for branches between joints

    labeled_segments = np.zeros_like(pruned_skeleton,dtype=int)
    assign_id = 1 
    visited_pixels = set()
    
    # If a specific starting pixel (northern most jointpoint) is provided, trace from it first
    if starting_pixel:
        labeled_segments = trace_branch_and_label(starting_pixel, starting_pixel, joint_coords, pruned_skeleton, labeled_segments, assign_id, visited_pixels)
        assign_id += 1
        # Remove starting pixel from joint coordinates to avoid re-tracing
        joint_coords = [coord for coord in joint_coords if tuple(coord) != tuple(starting_pixel)]


    # Iterate over joints to trace branches between joints
    for joint in joint_coords:
        neighbors = [(joint[0] + i, joint[1] + j)
                     for i in [-1, 0, 1] for j in [-1, 0, 1]
                     if not (i == 0 and j == 0)]
        
        branch_lengths = []
        branch_ids = []
        initial_joint = (joint[0],joint[1])
        
        
        for neighbor in neighbors:
            if pruned_skeleton[neighbor] == 1 and neighbor not in visited_pixels:
                # Trace the branch starting from this neighbor
                current_label_segments = np.copy(labeled_segments)
                current_label_segments, visited_pixels = trace_branch_and_label(neighbor, initial_joint, joint_coords, pruned_skeleton, current_label_segments, assign_id, visited_pixels)
                
                # Calculate branch length (for assignment of the joint pixel)
                branch_length = np.sum(current_label_segments == assign_id)
                if branch_length > 1:  # Proceed only if branch contains more than one pixel
                    branch_lengths.append(branch_length)
                    branch_ids.append(assign_id)
                    # Assign unique ID for each traced branch
                    labeled_segments = np.where(current_label_segments == assign_id, assign_id, labeled_segments)
                    assign_id += 1

            # Find longest branch connected to the joint and assign its ID to the joint pixel
            if branch_lengths:
                longest_branch_id = branch_ids[np.argmax(branch_lengths)]
                labeled_segments[joint[0], joint[1]] = longest_branch_id # Set joint ID to the ID of the longest connected branch
                visited_pixels.add((joint[0], joint[1]))


    # Assign ID to redundant joint pixels (branch len == 1)
    joint_coords_redundant = np.column_stack(np.where((pruned_skeleton == 1) & (labeled_segments == 0)))

    if len(joint_coords_redundant) > 0:
        labeled_coords = np.column_stack(np.where(labeled_segments > 0))
        branch_ids = labeled_segments[labeled_coords[:, 0], labeled_coords[:, 1]]
        branch_tree = cKDTree(labeled_coords)
        updated_labeled_segments = labeled_segments.copy()

        for joint in joint_coords_redundant:
            distance, nearest_index = branch_tree.query(joint)
            nearest_branch_pixel = labeled_coords[nearest_index]
            nearest_branch_id = labeled_segments[tuple(nearest_branch_pixel)]
            updated_labeled_segments[tuple(joint)] = nearest_branch_id
    else:
        updated_labeled_segments = labeled_segments.copy()


    return updated_labeled_segments


def order_raster_from_endpoint(subskel,branch):
    # - Input is a skeleton with exactly one branch, with value =1 on the branch and zeros elsewhere
    # - Traces the branch from the first identified endpoint to the last
    # - returns ordered pixel coordinates in two arrays (x,y)


    endpoints = find_endpoints(subskel)

    if len(np.where(endpoints==1)[0]) > 2:
        print('Caution! Branch No. '+str(branch)+', greater than 2 endpoints.')


    # Edge case when less than 2 endpoints are found
    if len(np.where(endpoints==1)[0]) < 2:
        print('Branch No. '+str(branch)+', less than 2! Skeletonize...')
        subskel = skeletonize(subskel)
        endpoints = find_endpoints(subskel)

        if len(np.where(endpoints==1)[0]) < 2:
            print('Branch No. '+str(branch)+', still less than two endpoints after skeletonization...try alternate route.')
        
            joints, joints_removed = find_joints(subskel)
            endpoint_candidates = find_endpoints(joints_removed)
  
            endpoint_candidates = np.where(endpoint_candidates & endpoints, False, endpoint_candidates)


            # Select pixels with non-direct neighbors
            indirect_neighbors = find_indirect_neighbors(joints_removed)
            endpoint_candidates = np.where(endpoint_candidates & indirect_neighbors, True, False)

            # Choose one! remove the other
            ep2y,ep2x = np.where(endpoint_candidates == 1)

            # Set second endpoint 
            endpoints[ep2y[0]][ep2x[0]] = 1

            # Remove other endpoint candidate from skeleton
            subskel[ep2y[1]][ep2x[1]] = 0

            plt.figure()
            plt.imshow(subskel)
            y,x = np.where(endpoint_candidates == 1)
            plt.scatter(x, y, color='red', s=0.5)
            y, x = np.where(endpoints == 1)
            plt.scatter(x, y, color='blue', s=0.5) 
            plt.xlim(1620,1630)
            plt.ylim(3150,3200)
            plt.show()

            print('Number of endpoints: ',len(np.where(endpoints == 1)[0]))


    firstpoint = (np.where(endpoints)[0][0], np.where(endpoints)[1][0]) # y, x
    lastpoint = [np.where(endpoints)[0][1], np.where(endpoints)[1][1]] # y, x

    
    # Now, from this point, trace neighbors and return ordered pixel_y and pixel_x list
    ordered_branch_pixels = trace_branch_from_endpoint(firstpoint,[lastpoint],subskel)
    
    # Append last point
    ordered_branch_pixels.append((np.where(endpoints)[0][1], np.where(endpoints)[1][1]))
    
    pixel_coords_x = [x for y,x in ordered_branch_pixels]
    pixel_coords_y = [y for y,x in ordered_branch_pixels]

    return pixel_coords_x,pixel_coords_y

# ---------------------------- (PRIMARY FUNCTIONS) ----------------------------


def get_skeleton(mask_image,pixcdate,figdir,dilate_amount,gauss_amount,savePlot=False):
    # Function for extracting the initial skeleton from binary raster image
    # INPUTS:
    #       - mask_image: binary raster image of water (1s) and not water (0s)
    #       - pixcdate: string of date 'yyyymmdd' of the pixel cloud acquisition
    #       - figdir: figure directory
    #       - dilate_amount: parameter controlling pixel dilation (meters)
    #       - gauss_amount: parameter controlling gaussian filter kernal size (meters)
    # OUTPUTS:
    #       - skeleton: binary image of returned skeleton

    print('...........Running get skeleton algorithm...........')


    LFmonths = [1,2,3,4,5,11,12]
    HFmonths = [6,7,8,9,10]

    if int(pixcdate[4:6]) in HFmonths:
        season = 'HF'
    if int(pixcdate[4:6]) in LFmonths:
        season = 'LF'
    

    if season == 'LF':
        # Step 1: buffer water mask to combine small geometries (across dams, bridges, etc)
        dilated_mask = binary_dilation(mask_image, footprint=createKernel(dilate_amount)) # OG 10

        # Step 2: Select only the largest connected structure
        labeled_mask, num_features = label(dilated_mask, return_num=True)
        component_sizes = np.bincount(labeled_mask.ravel())

        # - Ignore background (label 0)
        component_sizes[0] = 0

        # - Find the label of the largest component and select only mask with this label
        largest_component_mask = (labeled_mask == component_sizes.argmax())

        # Step 3. Perform gaussian smoothing
        smoothed_mask = gaussian(largest_component_mask.astype(float), sigma=float(gauss_amount)) # OG 20
        smoothed_mask = smoothed_mask > 0.6  # Adjustable threshold

        # Step 3.2: Perform secondary dilation
        dilated_mask2 = binary_dilation(smoothed_mask, footprint=createKernel(dilate_amount)) # OG 10

        # Step 3.1: Select again the largest connected structure
        labeled_mask, num_features = label(dilated_mask2, return_num=True)
        component_sizes = np.bincount(labeled_mask.ravel())
        component_sizes[0] = 0 # Setting largest component (being the background) to zero... need to make this more robust in case background is not connected
        largest_component_mask_smoothed = (labeled_mask == component_sizes.argmax())


    else:
        dilated_mask = binary_dilation(mask_image, footprint=createKernel(dilate_amount)) # OG 10

        # Step 2: Select only the largest connected structure
        labeled_mask, num_features = label(dilated_mask, return_num=True)
        component_sizes = np.bincount(labeled_mask.ravel())

        # - Ignore background (label 0)
        component_sizes[0] = 0

        # - Find the label of the largest component and select only mask with this label
        largest_component_mask = (labeled_mask == component_sizes.argmax())

        # Step 3. Perform gaussian smoothing
        smoothed_mask = gaussian(largest_component_mask.astype(float), sigma=float(gauss_amount)) # OG 30
        smoothed_mask = smoothed_mask > 0.6  # Adjustable threshold (OG 0.8)

        # Step 3.2: Perform secondary dilation
        dilated_mask2 = binary_dilation(smoothed_mask, footprint=createKernel(dilate_amount)) # OG 20

        # Step 3.1: Select again the largest connected structure
        labeled_mask, num_features = label(dilated_mask2, return_num=True)
        component_sizes = np.bincount(labeled_mask.ravel())
        component_sizes[0] = 0 # Setting largest component (being the background) to zero... need to make this more robust in case background is not connected
        largest_component_mask_smoothed = (labeled_mask == component_sizes.argmax())



    # Step 4. Perform skeletonization
    # With gaussian smoothing:
    skeleton = skeletonize(largest_component_mask_smoothed)

    if savePlot == True:
        # Step 5. Find endpoints
        endpoints = find_endpoints(skeleton)

        # display results and save to figure
        fig, axes = plt.subplots(nrows=1, ncols=6, figsize=(16, 4), sharex=True, sharey=True)
        ax = axes.ravel()

        ax[0].imshow(mask_image, cmap=plt.cm.gray)
        ax[0].axis('off')
        ax[0].set_title('Water mask', fontsize=15)

        ax[1].imshow(dilated_mask, cmap=plt.cm.gray)
        ax[1].axis('off')
        ax[1].set_title('Dilation', fontsize=15)

        ax[2].imshow(largest_component_mask, cmap=plt.cm.gray)
        ax[2].axis('off')
        ax[2].set_title('Largest connection', fontsize=15)

        ax[3].imshow(smoothed_mask, cmap=plt.cm.gray)
        ax[3].axis('off')
        ax[3].set_title('Gaussian filter', fontsize=15)

        ax[4].imshow(largest_component_mask_smoothed, cmap=plt.cm.gray)
        ax[4].axis('off')
        ax[4].set_title('Secondary dilation', fontsize=15)

        ax[5].imshow(np.zeros_like(skeleton), cmap='gray')
        y, x = np.where(skeleton == 1)  # Get skeleton coordinates
        ax[5].scatter(x, y, color='red', s=0.05)  # Plot skeleton points as red dots
        ax[5].axis('off')
        ax[5].set_title('Skeleton', fontsize=15)

        fig.tight_layout()
        #fig.suptitle("Dilation: "+str(dilate_amount)+", Gauss: "+str(gauss_amount))
        plt.savefig(figdir+'/'+pixcdate+'gauss_dilate_'+str(gauss_amount)+str(dilate_amount)+'extractedSkeleton.png')
        plt.close()

    return skeleton


def get_pruned_skeleton(skeleton,mask_image,pixcdate,figdir,distance_threshold=30,max_hole_size=300,prunethresh=600):
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
    print('...........Running branch pruning algorithm...........')


    # display results
    fig, axes = plt.subplots(nrows=1, ncols=6, figsize=(16, 4), sharex=True, sharey=True)
    ax = axes.ravel()

    ax[0].imshow(np.zeros_like(skeleton), cmap='gray')
    y, x = np.where(skeleton == 1) 
    ax[0].scatter(x, y, color='red', s=0.05)
    ax[0].axis('off')
    ax[0].set_title('Skeleton', fontsize=15)

    # Connect close endpoints
    # pixels (1 pixel = 10 m)
    connected_skeleton = connect_close_endpoints(skeleton, distance_threshold)

    # Compute joints and endpoints for plotting
    endpoints = find_endpoints(connected_skeleton)
    joints, joints_removed_mask = find_joints(connected_skeleton)


    ax[1].imshow(np.zeros_like(skeleton), cmap='gray')
    y, x = np.where(connected_skeleton == 1) 
    ax[1].scatter(x, y, color='red', s=0.05)  
    y, x = np.where(endpoints == 1) 
    ax[1].scatter(x, y, color='yellow', s=0.5) 
    y, x = np.where(joints == 1) 
    ax[1].scatter(x, y, color='blue', s=0.5) 
    ax[1].axis('off')
    ax[1].set_title('Connect, E&J', fontsize=15)


    #Fill small holes
    filled_skeleton = fill_small_holes(connected_skeleton, max_hole_size) # hole size is radius in pixels

    # skeletonize again
    skeleton2 = skeletonize(filled_skeleton)

    # Compute joints and endpoints for plotting
    endpoints = find_endpoints(skeleton2)
    joints, joints_removed_mask = find_joints(skeleton2)


    ax[2].imshow(np.zeros_like(skeleton), cmap='gray')
    y, x = np.where(filled_skeleton == 1) 
    ax[2].scatter(x, y, color='red', s=0.05)  
    y, x = np.where(skeleton2 == 1) 
    ax[2].scatter(x, y, color='yellow', s=0.5) 
    ax[2].axis('off')
    ax[2].set_title('Fill holes', fontsize=15)

    ax[3].imshow(np.zeros_like(skeleton), cmap='gray')
    y, x = np.where(skeleton2 == 1) 
    ax[3].scatter(x, y, color='red', s=0.05)  
    y, x = np.where(endpoints == 1) 
    ax[3].scatter(x, y, color='yellow', s=0.5) 
    y, x = np.where(joints == 1) 
    ax[3].scatter(x, y, color='blue', s=0.5) 
    ax[3].axis('off')
    ax[3].set_title('E&J', fontsize=15)

    
    # -------------------- Recompute endpoints/jointpoints, then iterate through pruning algorithm until no changes occur!
    pruned_skeleton = prune_branches(skeleton2, prunethresh)
    pruned_skeleton = skeletonize(pruned_skeleton)
    endpoints_len_new = float('inf')
    i = 0
    while True:
        endpoints = find_endpoints(pruned_skeleton)
        y, x = np.where(endpoints == 1) 
        endpoints_len = len(x)

        if endpoints_len_new == endpoints_len:
            final_skeleton = skeletonize(pruned_skeleton) # get rid of areas with 2 pixels width
            break
        
        print('Running prune branches...')
        pruned_skeleton = prune_branches(pruned_skeleton, prunethresh) # APPLY pruning algorithm
        i +=1
        # Compute length of endpoints
        y, x = np.where(endpoints == 1) 
        endpoints_len_new = len(x)

    print('Pruning algorithm applied '+str(i)+' times!')


    ax[4].imshow(np.zeros_like(skeleton), cmap='gray')
    y, x = np.where(skeleton2 == 1)
    ax[4].scatter(x, y, color='red', s=0.05)
    y, x = np.where(final_skeleton == 1) 
    ax[4].scatter(x, y, color='yellow', s=0.5) 
    ax[4].axis('off')
    ax[4].set_title('Prune branches', fontsize=15)

    # Recompute joints and ends for plotting
    endpoints = find_endpoints(final_skeleton)
    joints, joints_removed_mask = find_joints(final_skeleton)

    ax[5].imshow(mask_image, cmap=plt.cm.gray)
    y, x = np.where(pruned_skeleton == 1)
    ax[5].scatter(x, y, color='red', s=0.05)
    y, x = np.where(endpoints == 1) 
    ax[5].scatter(x, y, color='yellow', s=0.5) 
    y, x = np.where(joints == 1) 
    ax[5].scatter(x, y, color='blue', s=0.5) 
    ax[5].axis('off')
    ax[5].set_title('Final skeleton', fontsize=15)

    fig.tight_layout()
    plt.savefig(figdir+'/'+pixcdate+'prunedSkeleton.png')
    #plt.close()


    return final_skeleton


def get_labeled_skeleton(final_skeleton,pixcdate,figdir,savePlot=False):

    print('...........Running skeleton labeling algorithm...........')

    # Now, set ID values to each channel....
    endpoints = find_endpoints(final_skeleton)
    joints, joints_removed_mask = find_joints(final_skeleton)
    # id = np.where(endpoints)[0] == np.min(np.where(endpoints)[0])
    # start_pixel = (np.where(endpoints)[0][id][0], np.where(endpoints)[1][id][0]) # (y,x), choose endpoint with the smallest y value (most north endpoint)

    joint_coords = np.where(joints == 1)  # Get skeleton coordinates
    joint_coords = np.dstack((joint_coords[0],joint_coords[1]))
    joint_coords = joint_coords.reshape(-1, 2)

    # Assign labels
    labeled_skeleton = assign_unique_ids_to_branches(final_skeleton,joint_coords,None)

    if savePlot == True:
        # Get coordinates and labels of all labeled pixels
        y_coords, x_coords = np.where(labeled_skeleton > 0)
        labels = labeled_skeleton[y_coords, x_coords]

        # display results
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(4, 4), sharex=True, sharey=True)

        ax = axes.ravel()

        ax[0].imshow(np.zeros_like(final_skeleton), cmap='gray')
        y, x = np.where(final_skeleton == 1)  # Get skeleton coordinates
        ax[0].scatter(x, y, color='red', s=0.05)  # Plot skeleton points as red dots
        y, x = np.where(endpoints == 1) 
        ax[0].scatter(x, y, color='yellow', s=0.5) 
        y, x = np.where(joints == 1) 
        ax[0].scatter(x, y, color='blue', s=0.5) 
        ax[0].axis('off')
        ax[0].set_title('E&J', fontsize=15)


        # View label
        ax[1].imshow(np.zeros_like(final_skeleton), cmap='gray')
        ax[1].scatter(x_coords, y_coords, c=labels, cmap='tab20', s=1)

        # Annotate each branch with its ID near its center
        unique_labels = np.unique(labels)
        for branch_id in unique_labels:
            branch_yx = np.column_stack(np.where(labeled_skeleton == branch_id))
            if len(branch_yx) > 0:
                # Find approximate center by using the median position
                center_y, center_x = np.median(branch_yx, axis=0).astype(int)
                ax[1].text(center_x+300, center_y, str(branch_id), color="white", fontsize=8, ha='center')

        ax[1].axis('off')
        y, x = np.where((final_skeleton == 1) & (labeled_skeleton == 0)) 
        ax[1].scatter(x, y, color='yellow', s=1) # These are skeleton points without a label


        ax[1].set_title('Labeled branch', fontsize=15)

        fig.tight_layout()
        plt.savefig(figdir+'/'+pixcdate+'labeledSkeleton.png')
        plt.close()


    return labeled_skeleton


def extract_cl_from_skeleton(labeled_skeleton,water_mask_tiff):

    print('...........Running centerline extraction from skeleton algorithm...........')

    line_data = []
    for branch in np.unique(labeled_skeleton):
        if branch == 0:
            continue

        # extract subskeleton
        subskel = np.zeros_like(labeled_skeleton)
        subskel[labeled_skeleton == branch] =  1

        # order raster from endpoint
        pixel_coords_x,pixel_coords_y = order_raster_from_endpoint(subskel,branch)

        # Convert ordered pixel coordinates to lat lon
        line = pixel_coordinates2latlonline(water_mask_tiff,pixel_coords_y,pixel_coords_x,branch,reverse=True)

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

def plot_merged(cl_merged, maskdate,figdir):

    # Ensure GeoDataFrame is in a projected CRS for accurate placement
    if not cl_merged.crs.is_projected:
        print("Reprojecting to UTM for better plotting.")
        cl_merged = cl_merged.to_crs('EPSG:3857')  # Example of a projected CRS

    # Create a color map for branch IDs
    unique_ids = cl_merged['branch_id'].unique()
    num_branches = len(unique_ids)
    # Generate a repeating discrete color map using tab20 or a larger palette
    colors = plt.cm.tab20(np.linspace(0, 1, 20))  # tab20 has 20 colors
    color_map = {branch_id: colors[i % 20] for i, branch_id in enumerate(unique_ids)}

    # Plot the linestrings
    fig, ax = plt.subplots(figsize=(3, 6))
    for _, row in cl_merged.iterrows():
        color = color_map[row['branch_id']]
        ax.plot(*row.geometry.xy, color=color, linewidth=2, label=f"Branch {row['branch_id']}")

    # Remove duplicate labels in the legend
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #ax.legend(by_label.values(), by_label.keys(), loc='upper left', bbox_to_anchor=(1, 1))

    # Add text labels for each branch ID
    for _, row in cl_merged.iterrows():
        x, y = row.geometry.xy[0][len(row.geometry.xy[0]) // 2], row.geometry.xy[1][len(row.geometry.xy[1]) // 2]
        ax.text(x, y, str(row['branch_id']), fontsize=5, ha='center', va='center', color='black')

    # Set plot properties
    ax.set_title("Generated "+str(maskdate))
    ax.axis("equal")
    # Remove x and y tick marks
    ax.set_xticks([])
    ax.set_yticks([])

    plt.tight_layout()
    plt.savefig(figdir+'/'+maskdate+'_generated_centerlines.png')
    plt.show()


# ---------------------------- (SECONDARY FUNCTIONS: Used for plotting/merging centerline portions to show slope in connected channels) ----------------------------


def fill_cl_gaps(sel_cl,hemi,distance_threshold = 1000):
    
    # Input is a geodataframe of reach centerlines, each with a single linestring geometry
    # This function creates a new linestring geometry to bridge all gaps greater than 1 km

    centroid = sel_cl.geometry.unary_union.centroid
    utm_zone = int((centroid.x + 180) // 6) + 1  # Calculate UTM zone
    if hemi == 'north':
        utm_crs = f'EPSG:{32600 + utm_zone}'
    if hemi == 'south':
        utm_crs = f'EPSG:{32700 + utm_zone}'

    sel_cl_utm = sel_cl.to_crs(utm_crs)

    endpoints = []
    for idx, line in sel_cl_utm.iterrows():
        if isinstance(line.geometry, LineString):
            start_point = Point(line.geometry.coords[0])
            end_point = Point(line.geometry.coords[-1])

            endpoints.append({"branch_id": line.branch_id, "geometry": start_point, "type": "start"})
            endpoints.append({"branch_id": line.branch_id, "geometry": end_point, "type": "end"})

    endpoints_gdf = gpd.GeoDataFrame(endpoints, crs=sel_cl_utm.crs)

    # Extract coordinates of endpoints
    coords = np.array([(point.x, point.y) for point in endpoints_gdf.geometry])
    tree = cKDTree(coords)

    # Query nearest neighbors
    pairs = tree.query_pairs(distance_threshold)

    # Create connections for closest endpoint pairs
    connections = []
    for idx1, idx2 in pairs:
        point1 = endpoints_gdf.iloc[idx1].geometry
        point2 = endpoints_gdf.iloc[idx2].geometry
        connections.append(LineString([point1, point2]))

    # Create a GeoDataFrame for connections
    connections_gdf = gpd.GeoDataFrame(geometry=connections, crs=sel_cl_utm.crs)

    # Combine original geometries and connections
    merged_gdf = gpd.GeoDataFrame(
        pd.concat([sel_cl_utm, connections_gdf], ignore_index=True), 
        crs=sel_cl_utm.crs
    )
    merged_gdf = merged_gdf.to_crs('EPSG:4326')

    # check that merged_gdf can be merged to a single linestring
    multi_line = MultiLineString(list(merged_gdf.geometry))

    # # Merge only connected parts
    merged_parts = linemerge(multi_line)

    # Check to see if all lines in multilinestring are significant
    newline = []
    minlength = 0.001

    if isinstance(merged_parts,LineString):
            return merged_parts
    else:

        for line in list(merged_parts.geoms):
            # check length of line, ignore if less tham 0.001 deg
            if line.length > minlength:
                newline.append(line)

        merged_parts = MultiLineString(newline)

        return merged_parts
    


def merge_multiline_with_gaps(multi_line):
    """
    For use when 'fill_cl_gaps' misses a section
    This function created with help of ChatGPT. prompt: 
    I need to be able to input a multilinestring that contains 2 or more lines, 
    identify the gaps between these lines (where the lines are closest together), 
    and bridge them to produce one single linestring
    """
    if len(multi_line.geoms) == 1:
        # If there's only one LineString, return it directly
        return multi_line.geoms[0]
    
    # List to store the connecting LineStrings
    connectors = []
    remaining_lines = list(multi_line.geoms)

    # Keep bridging gaps until all LineStrings are connected
    while len(remaining_lines) > 1:
        # Initialize variables to store the closest pair and the minimum distance
        min_distance = float("inf")
        closest_pair = None

        # Compare each pair of LineStrings to find the closest points
        for i, line1 in enumerate(remaining_lines):
            for j, line2 in enumerate(remaining_lines):
                if i >= j:  # Skip duplicate and self-pairs
                    continue
                point1, point2 = nearest_points(line1, line2)
                distance = point1.distance(point2)
                if distance < min_distance:
                    min_distance = distance
                    closest_pair = (i, j, point1, point2)

        # Extract the closest pair and their points
        i, j, point1, point2 = closest_pair
        connector = LineString([point1, point2])
        connectors.append(connector)

        # Merge the two closest lines and replace them in the remaining list
        merged_line = linemerge(MultiLineString([remaining_lines[i], connector, remaining_lines[j]]))
        remaining_lines = [line for k, line in enumerate(remaining_lines) if k not in (i, j)]
        remaining_lines.append(merged_line)

    # Return the final connected LineString
    return remaining_lines[0]


def add_line_geom(start, end):
    # Function to add a straight line between two coordinate points
    y0, x0 = start
    y1, x1 = end
    y_values = np.array([y0,y1])
    x_values = np.array([x0,x1])

    points = list(zip(y_values, x_values))
    line = LineString(points)
    return line
