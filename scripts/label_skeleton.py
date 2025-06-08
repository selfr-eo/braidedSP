import glob
import sys
import os
from tqdm import tqdm

import geopandas as gpd

from braided_rivers.braided import Mask
import skeletonize_func as skel

# ------------------------------ Paths ------------------------------

# read shapefile for starting line
river_bounds = gpd.read_file(r"C:\Users\cwch\Projects\eo4flood_local\Niger_AOI_DHI1D.shp")
river_crs = river_bounds.crs
river_bounds = river_bounds.loc[0, 'geometry']

# read the starting line for labeling
starting_line = gpd.read_file(r"\\dkcph1-nas02\jupyterhub-exchange\connorchewning\data\eo4flood\Niger_Hydrodynamic\niger_river_centerline_processing\niger_starting_line.shp").to_crs(river_crs)
main_centerline = gpd.read_file(r"C:\Users\cwch\Projects\eo4flood_local\niger_2d_area_SWORD.shp")
main_centerline = skel.trim_cl_to_river_bounds(main_centerline.to_crs(river_crs), river_bounds)


# read in originally processed centerline
cl = gpd.read_file(r"C:\Users\cwch\Projects\eo4flood_local\Niger_2D_cl_merged.geojson")
cl = skel.trim_cl_to_river_bounds(cl.to_crs(river_crs), river_bounds)

# joint at the joints
cl = skel.join_cl_at_joints(cl, starting_line, main_centerline, search_dist=40)
cl.to_file(r"C:\Users\cwch\Projects\eo4flood_local\Niger_2D_cl_merged_adjusted.geojson")