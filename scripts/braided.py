# ------------------------------ Imports ------------------------------
from dataclasses import dataclass

import os


# import numpy as np
# import pandas as pd
# import geopandas as gpd
# import matplotlib.pyplot as plt

# # Computing slopes
# from sklearn.linear_model import LinearRegression
# from shapely.geometry import LineString
# from skimage.filters import gaussian

# # import custom functions
# import procBraided as pc
import skeletonize_func as skel


@dataclass
class Mask:
    path: str
    odir: str
    identifier: str

    def __post_init__(self):
        # make sure file exists
        if not os.path.exists(self.path):
            raise FileNotFoundError(
                "A mask can only be created on existing raster files."
            )

        # set up export dirs
        if not os.path.exists(self.odir):
            os.makedirs(self.odir)

        # set up fig directory
        self.figdir = os.path.join(self.odir, "figs")
        if not os.path.exists(self.figdir):
            os.makedirs(self.figdir)

        # read in the file to memory
        self.array = skel.get_watermask(self.path)

    def process(
        self,
        dilate=10,
        gauss=0.5,
        fill_hole_size=100,
        hemi="north",
        season="HF",
        second_dilation=False,
        distance_threshold=5,
        prune_thresh=600,
        river_bounds=None,
        save_progress=False,
    ):
        
        # smooth the mask
        smoothed_mask = skel.smooth_mask(
            self.array,
            self.identifier,
            self.figdir,
            dilate,
            gauss,
            max_hole_size = fill_hole_size,
            season=season,
            second_dilation=second_dilation,
            savePlot=False,
        )
        if save_progress:
            skel.export_progress(smoothed_mask, self.path, os.path.join(self.figdir, f"smoothed_mask_{self.identifier}.tif"))


        # gets pruned skeleton from the mask
        skeleton = skel.get_skeleton(
            smoothed_mask,
            self.array,
            self.identifier,
            self.figdir,
            distance_threshold=distance_threshold, # distance thresh for connecting pixels
            prunethresh=prune_thresh,
        )  
        if save_progress:
            skel.export_progress(skeleton, self.path, os.path.join(self.figdir, f"{self.identifier}_skeleton.tif"))


        # labels the reaches from top to bottom (improvements can be made here by letting user provide a starting point)
        labeled_skeleton = skel.get_labeled_skeleton(
            skeleton, self.identifier, self.figdir, savePlot=False
        )
        if save_progress:
            skel.export_progress(labeled_skeleton, self.path, os.path.join(self.figdir, f"{self.identifier}_labeled_skeleton.tif"))


        # take the results of processing the skeleton and generate centerlines
        cl = skel.extract_cl_from_skeleton(labeled_skeleton, self.path)

        # trim the centerlines based on a starting point for the river system
        if river_bounds:
            cl = skel.trim_cl_to_river_bounds(cl, river_bounds)

        # merge the centerline joints into intersections

        # Option to merge short centerlines, but will remove relevance of joints
        # cl = skel.merge_short_centerlines(cl, hemi) # options to merge short centerlines into a longer stretch
        self.cl = cl

    def export(self, path):
        self.cl.to_file(path)
