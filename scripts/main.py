# ------------------------------ Imports ------------------------------
import os
import glob
import sys

from tqdm import tqdm

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

# Computing slopes
from sklearn.linear_model import LinearRegression
from shapely.geometry import LineString
from skimage.filters import gaussian

# import custom functions
import procBraided as pc
import skeletonize_func as skel


# ------------------------------ Paths ------------------------------
# hard code paths for now / later add in as system or function arguments ?
maskdir = r"../data/watermasks"  # or later just have use provide paths to mask
mask_paths = glob.glob(os.path.join(maskdir, "*.tif"))
if len(mask_paths) == 0:
    raise ValueError("No mask files found in mask directory")

# set up export dirs
odir = r"../data/odir/"
if not os.path.exists(odir):
    os.makedirs(odir)

figdir = r"../data/odir/figs/"
if not os.path.exists(figdir):
    os.makedirs(figdir)


# ------------------------------ Parameters ------------------------------
hemi = "north"
dilate = 1 #pixels S2 res, so in 10 meters
gauss = 5
fill_hole_size = 0 # 0 means no hole filling


# ------------------------------ Main ------------------------------
# generate the centerlines for each mask file
for file in tqdm(mask_paths, desc="Processing masks"):

    ##### hardcoded date format to remove
    maskdate = os.path.basename(file)[24:28]
    pixcdate = (
        "2021" + maskdate[0:2] + "01"
    )  # format expected by the function is YYYYMMDD

    # reads the watermask geotiff
    water_mask = skel.get_watermask(file)


    # gets initial skeleton
    skeleton = skel.get_skeleton(
        water_mask,
        pixcdate,
        figdir,
        dilate,
        gauss,
        season="HF",
        second_dilation=False,
        savePlot=False,
    )

    # copy raster for exporting process updates
    skel.export_progress(skeleton, file, os.path.join(odir, maskdate + f"_dilate_{dilate}.tif"))
    #skel.export_progress(skeleton, file, os.path.join(odir, maskdate + f"_gauss_{gauss}.tif"))

    # prunes skeleton to ensure each reach has an inflow and outflow, also closes holes
    pruned_skeleton = skel.get_pruned_skeleton(
        skeleton,
        water_mask,
        pixcdate,
        figdir,
        distance_threshold=5,
        max_hole_size=fill_hole_size,
        prunethresh=600,
    )  # distance thresh for connecting pixels

    # labels the reaches from top to bottom (improvements can be made here by letting user provide a starting point)
    labeled_skeleton = skel.get_labeled_skeleton(
        pruned_skeleton, pixcdate, figdir, savePlot=True
    )

    # take the results of processing the skeleton and generate centerlines
    cl = skel.extract_cl_from_skeleton(labeled_skeleton, file)
    # cl = skel.merge_short_centerlines(cl, hemi) # options to merge short centerlines into a longer stretch

    # skel.plot_merged(cl, maskdate, figdir)

    # print("SUCCESSFUL SKELETON! Saving files...")
    # save processed shape 
    cl.to_file(
        odir + "/" + maskdate + "_fs_" + str(fill_hole_size) + f"_generated_cl_{dilate}_{gauss}_{fill_hole_size}.geojson"
    )

    sys.exit()