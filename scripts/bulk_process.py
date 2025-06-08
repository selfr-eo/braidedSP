import glob
import sys
import os
from tqdm import tqdm

import geopandas as gpd

from braided_rivers.braided import Mask

"""Script assumes that masksa are part of single river system and that parameters for processing are provided within the shape files
"""


# ------------------------------ Paths ------------------------------
maskdir = r"\\dkcph1-nas02\jupyterhub-exchange\connorchewning\data\eo4flood\Niger_Hydrodynamic\niger_river_centerline_processing\watermasks"  # or later just have use provide paths to mask
odir = r"\\dkcph1-nas02\jupyterhub-exchange\connorchewning\data\eo4flood\Niger_Hydrodynamic\niger_river_centerline_processing\odir"
river_chunk_path = r"\\dkcph1-nas02\jupyterhub-exchange\connorchewning\data\eo4flood\Niger_Hydrodynamic\niger_river_centerline_processing\niger_aoi_processing_sections.shp"

# ------------------------------ Parameters ------------------------------
chunck_id_col = 'FID'
area_id = "Niger_2D"
watermask_crs = "EPSG:4326"

# ------------------------------ Main ------------------------------
# generate the centerlines for each mask file
river_chunks = gpd.read_file(river_chunk_path).to_crs(watermask_crs)

# cycle through
for i in tqdm(river_chunks.index, desc="Processing masks"):

    # see if we need to process this chunk
    process = bool(river_chunks.loc[i, 'process'])
    if not process:
        break

    # Pull out chunk infor
    chunk_id = river_chunks.loc[i, chunck_id_col]
    river_bounds = river_chunks.loc[i, 'geometry']
    file = os.path.join(maskdir, f"{i}.tif")

    # create mask object
    params={
        'dilate': int(river_chunks.loc[i, 'dilate']),
        'gauss' : river_chunks.loc[i, 'gauss'],
        'fill_hole_size' : river_chunks.loc[i, 'fill_size']
    }
    identifier = f"{area_id}_chunk_{chunk_id}_dil{params['dilate']}_gauss{params['gauss']}_fill_{params['fill_hole_size']}"
    mask = Mask(file, odir, identifier)

    # process the mask to centerline skeleton
    mask.process(
        dilate = params['dilate'], #pixels S2 res, so in 10 meters
        gauss = params['gauss'], # standard deviation of gaussian filter
        fill_hole_size = params['fill_hole_size'], # 0 means no hole filling
        river_bounds=river_bounds, # a shapely object in the same crs as the tiff files
        save_progress=True
    )
    mask.export(os.path.join(odir, f"{identifier}_cl_labeled.geojson"))

    
    # ## TODO Move to mask functions in main package
    # # adjust the trimmed centerline to get conencted and labeled joints
    # import skeletonize_func as skel

    # # read in originally processed centerline
    # cl = gpd.read_file(rf"C:\Users\cwch\Tools\braided_rivers\data\odir\{identifier}_cl_labeled.geojson")
    # cl = skel.trim_cl_to_river_bounds(cl, river_bounds)

    # starting_line = gpd.read_file(r"C:\Users\cwch\Tools\braided_rivers\data\starting_line.shp")
    # main_centerline = gpd.read_file(r"C:\Users\cwch\Tools\braided_rivers\data\main_centerline.shp")
    # main_centerline = skel.trim_cl_to_river_bounds(main_centerline, river_bounds)

    # cl = skel.join_cl_at_joints(cl, starting_line, main_centerline, search_dist=40)
    # cl.to_file(rf"C:\Users\cwch\Tools\braided_rivers\data\odir\{identifier}_cl_updated.geojson")
    