import glob
import sys
import os
from tqdm import tqdm

import geopandas as gpd

from braidedSP.mask import Mask



# ------------------------------ Paths ------------------------------
# hard code paths for now / later add in as system or function arguments ?
maskdir = r"\\dkcph1-nas02\jupyterhub-exchange\connorchewning\data\eo4flood\Niger_Hydrodynamic\niger_river_centerline_processing\watermasks"  # or later just have use provide paths to mask
mask_paths = glob.glob(os.path.join(maskdir, "*.tif"))
if len(mask_paths) == 0:
    raise ValueError("No mask files found in mask directory")

# set up export dirs
odir = r"\\dkcph1-nas02\jupyterhub-exchange\connorchewning\data\eo4flood\Niger_Hydrodynamic\niger_river_centerline_processing\odir"

# read shapefile for starting line
river_bounds = gpd.read_file(r"C:\Users\cwch\CNR\EO4FLOOD-grp - Niger\08_Models\Hydraulic model\DHI\Niger_AOI_DHI1D.shp")
river_bounds = river_bounds.loc[0, 'geometry']

# ------------------------------ Main ------------------------------
# generate the centerlines for each mask file
for file in tqdm(mask_paths, desc="Processing masks"):

    # extract infor based on known file structure
    identifier = file.split('S2_Niger_WaterMask_')[-1].split('.tif')[0]

    params={
        'dilate': 2,
        'gauss' : 0.5,
        'fill_hole_size' : 1000
    }
    identifier = f"{identifier}_dil{params['dilate']}_gauss{params['gauss']}_fill_{params['fill_hole_size']}"

    mask = Mask(file, odir, identifier)    
    mask.process(
        dilate = params['dilate'], #pixels S2 res, so in 10 meters
        gauss = params['gauss'], # standard deviation of gaussian filter
        fill_hole_size = params['fill_hole_size'], # 0 means no hole filling
        river_bounds=river_bounds, # a shapely object in the same crs as the tiff files
        save_progress=True
    )
    mask.export(os.path.join(odir, f"{identifier}_cl_labeled.geojson"))

    ## TODO Move to mask functions in main package
    # adjust the trimmed centerline to get conencted and labeled joints
    import skeletonize_func as skel

    # read in originally processed centerline
    cl = gpd.read_file(rf"C:\Users\cwch\Tools\braided_rivers\data\odir\{identifier}_cl_labeled.geojson")
    cl = skel.trim_cl_to_river_bounds(cl, river_bounds)

    starting_line = gpd.read_file(r"C:\Users\cwch\Tools\braided_rivers\data\starting_line.shp")
    main_centerline = gpd.read_file(r"C:\Users\cwch\Tools\braided_rivers\data\main_centerline.shp")
    main_centerline = skel.trim_cl_to_river_bounds(main_centerline, river_bounds)

    cl = skel.join_cl_at_joints(cl, starting_line, main_centerline, search_dist=40)
    cl.to_file(rf"C:\Users\cwch\Tools\braided_rivers\data\odir\{identifier}_cl_updated.geojson")