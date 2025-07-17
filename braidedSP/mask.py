# Typing imports
from dataclasses import dataclass

# general imports
import os
from tqdm import tqdm

# functional imports
from datetime import datetime
import numpy as np
import geopandas as gpd
import rasterio
from skimage.filters import gaussian
from skimage.morphology import skeletonize, binary_dilation, binary_erosion, remove_small_holes
#from braidedSP import tools
import tools

@dataclass
class Mask:

    river_name: str
    date: datetime
    path: str
    odir: str

    def __post_init__(self):
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

        # read in the file to memory
        self.array, self.ref_data, self.crs, self.profile = _get_watermask(self.path)

        self.smoothed = None
        self.skeleton = None # this must be processed before it is set

    def get_date(self):
        return self.date.strftime("%Y-%m-%d")

    def process(
        self,
        dilate=10,
        gauss=0.5,
        fill_hole_size=100,
        second_dilation=False,
        distance_threshold=5,
        prune_thresh=600,
        show_progress=False,
        save_progress=False
    ):
        
        # smooth the mask
        self.smoothed = _smooth(
            self.array,
            dilate,
            gauss,
            max_hole_size = fill_hole_size,
            second_dilation=second_dilation,
            show_progress=show_progress
        )
        if save_progress:
            self.export_progress(self.smoothed, f"{self.river_name}_{self.get_date()}_smoothed_mask.tif")


        # gets pruned skeleton from the mask
        skeleton = _skeletonize(
            self.smoothed,
            distance_threshold=distance_threshold, # distance thresh for connecting pixels
            prunethresh=prune_thresh,
            show_progress=show_progress
        )
        if save_progress:
            self.export_progress(skeleton, f"{self.river_name}_{self.get_date()}_skeleton.tif")

        # labels the reaches from top to bottom (improvements can be made here by letting user provide a starting point)
        labeled_skeleton = _label(skeleton)
        if save_progress:
            self.export_progress(labeled_skeleton, f"{self.river_name}_{self.get_date()}_labeled_skeleton.tif")

        # set the skeleton attribute
        self.skeleton = labeled_skeleton



    def vectorize(self, show_progress=False):

        labeled_skeleton = self.skeleton.copy()

        line_data = []
        skipped_branches = []
        for branch in tqdm(np.unique(labeled_skeleton), desc='Vectorizing branches',  leave=False, disable=not show_progress):
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
            line = tools.pixel_coordinates2latlonline(self.path, pixel_coords_y,pixel_coords_x,branch,reverse=True)

            # append gdf line to overall gdf with branch id
            line_data.append({'branch_id': branch, 'geometry': line})

        # go through skipped branches, and append to nearest branch
        for branch in tqdm(skipped_branches, desc='Vectorizing skipped branches',  leave=False, disable=not show_progress):

            # extract subskeleton
            subskel = np.zeros_like(labeled_skeleton)
            subskel[labeled_skeleton == branch] =  1

            # locate nearest pixel from another branch
            nearest_branch = _find_nearest_branch(subskel,labeled_skeleton)

            # merge pixels
            labeled_skeleton[labeled_skeleton == branch] = nearest_branch

            # rerun order raster from endpoint and convert to lat lon for the nearest branch with new appended pixel
            subskel = np.zeros_like(labeled_skeleton)
            subskel[labeled_skeleton == nearest_branch] =  1
            pixel_coords_x,pixel_coords_y = tools.order_raster_from_endpoint(subskel,nearest_branch)

            # Convert ordered pixel coordinates to lat lon
            line = tools.pixel_coordinates2latlonline(self.path, pixel_coords_y, pixel_coords_x,nearest_branch, reverse=True)

            # REMOVE old line with branch_id = nearest_branch from line_data
            line_data = [line for line in line_data if line['branch_id'] != nearest_branch]

            # append gdf line to overall gdf with branch id
            line_data.append({'branch_id': nearest_branch, 'geometry': line})

        cl_gdf_labeled = gpd.GeoDataFrame(line_data, crs=self.crs)

        return cl_gdf_labeled


    def export_progress(self, array, file_name):
        export_file = os.path.join(self.figdir, file_name)
        with rasterio.open(export_file, "w", **self.profile) as dst:

            ref_data = self.ref_data
            ref_data[0] = array
            dst.write(ref_data)






# ---------------------------- (Helper FUNCTIONS) ----------------------------
def _get_watermask(path):

    # Extract water mask data from provided geotiff using rasterio
    with rasterio.open(path) as src:
        water_mask_data = src.read(1)  # Read the first band (binary mask)
        ref_data = src.read()
        crs = src.crs
        profile = src.profile

    return water_mask_data, ref_data, crs, profile



def _smooth(mask_image, dilate_amount, gauss_amount, gauss_thresh=0.6, max_hole_size=300, second_dilation=True, show_progress=False):
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
    bar = tqdm(total=7, desc="Smoothing mask", leave=False, disable=not show_progress)

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




def _skeletonize(smoothed_mask, distance_threshold=30, prunethresh=600, show_progress=False):

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
    bar = tqdm(total=4, desc="Pruning", leave=False, disable=not show_progress)

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


def _label(final_skeleton):
    
    joints, joints_removed_mask = tools.find_joints(final_skeleton)

    joint_coords = np.where(joints == 1)  # Get skeleton coordinates
    joint_coords = np.dstack((joint_coords[0],joint_coords[1]))
    joint_coords = joint_coords.reshape(-1, 2)

    # Assign labels
    labeled_skeleton = tools.assign_unique_ids_to_branches(final_skeleton, joint_coords, None)

    return labeled_skeleton




def _find_nearest_branch(subskel,labeled_skeleton):

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
