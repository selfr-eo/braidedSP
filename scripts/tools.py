# Import libraries
import os
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, LineString, MultiLineString
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


# ---------------------------- (HELPER FUNCTIONS) ----------------------------
# (Functions in order of appearance)

def createKernel(radius):
    # create circular kernel for gaussian smoothing
    kernel = np.zeros((2*radius+1, 2*radius+1))
    y,x = np.ogrid[-radius:radius+1, -radius:radius+1]
    mask = x**2 + y**2 <= radius**2
    kernel[mask] = 1
    return kernel

def find_largest_component(mask):

    # take a mask, label the connected components and keep only the largest one
    labeled_mask, num_features = label(mask, return_num=True)
    component_sizes = np.bincount(labeled_mask.ravel())
    component_sizes[0] = 0 # Setting largest component (being the background) to zero... need to make this more robust in case background is not connected
    largest_component = (labeled_mask == component_sizes.argmax())

    return largest_component


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


def trace_branch_from_endpoint(endpoint, joint_coords, pruned_skeleton):
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
                # print('No pixels in branch!')
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
    for joint in tqdm(joint_coords, desc="Labeling", leave=False):
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


def order_raster_from_endpoint(subskel, branch):
    # - Input is a skeleton with exactly one branch, with value =1 on the branch and zeros elsewhere
    # - Traces the branch from the first identified endpoint to the last
    # - returns ordered pixel coordinates in two arrays (x,y)


    endpoints = find_endpoints(subskel)

    # if len(np.where(endpoints==1)[0]) > 2:
    #     print('Caution! Branch No. '+str(branch)+', greater than 2 endpoints.')


    # Edge case when less than 2 endpoints are found
    if len(np.where(endpoints==1)[0]) < 2:
        # print('Branch No. '+str(branch)+', less than 2! Skeletonize...')
        subskel = skeletonize(subskel)
        endpoints = find_endpoints(subskel)

        if len(np.where(endpoints==1)[0]) < 2:
            # print('Branch No. '+str(branch)+', still less than two endpoints after skeletonization...try alternate route.')
        
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
            #plt.show()
            plt.close()

            # print('Number of endpoints: ',len(np.where(endpoints == 1)[0]))


    firstpoint = (np.where(endpoints)[0][0], np.where(endpoints)[1][0]) # y, x
    lastpoint = [np.where(endpoints)[0][1], np.where(endpoints)[1][1]] # y, x

    
    # Now, from this point, trace neighbors and return ordered pixel_y and pixel_x list
    ordered_branch_pixels = trace_branch_from_endpoint(firstpoint,[lastpoint],subskel)
    
    # Append last point
    ordered_branch_pixels.append((np.where(endpoints)[0][1], np.where(endpoints)[1][1]))
    
    pixel_coords_x = [x for y,x in ordered_branch_pixels]
    pixel_coords_y = [y for y,x in ordered_branch_pixels]

    return pixel_coords_x,pixel_coords_y


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
