# Typing imports
import os
from dataclasses import dataclass
from datetime import datetime

# import project specific objects
from braidedSP.mask import Mask
from braidedSP.centerline import Centerline
from braidedSP.swot import SWOT

# general imports
from tqdm import tqdm
import pandas as pd
import geopandas as gpd


@dataclass
class River:

    name: str
    aoi: gpd.GeoDataFrame
    outdir: str

    def __post_init__(self):

        self.masks = []
        self.centerlines = []
        self.swot_obs = []

    def add_mask(self, mask_path:str, date:datetime):

        processed_mask = Mask(
            river_name = self.name,
            date = date,
            path = mask_path,
            odir = self.outdir,
        )

        self.masks.append(processed_mask)

    def process_masks(self, kwargs):

        # Progress bar settings
        show_progress = False
        if 'show_progress' in kwargs:
            show_progress = kwargs['show_progress']

        # see if we get multiple dictionaries as input
        if type(kwargs) is dict():
            kwargs = [kwargs] * len(self.masks)


        # cycle through and process masks
        for mask, mask_args in tqdm(zip(self.masks, kwargs), desc='Processing masks', leave=True, disable=not show_progress):
            mask.process(**mask_args)

    def generate_centerlines(self, kwargs):

        # Progress bar settings
        show_progress = False
        if 'show_progress' in kwargs:
            show_progress = kwargs['show_progress']

        # cycle through and process masks
        for mask in tqdm(self.masks, desc='Generating centerlines', leave=True, disable=not show_progress):

            # generate the centerline object and geometries
            gdf = mask.vectorize(**kwargs)
            cl = Centerline(
                river_name = mask.river_name,
                date = mask.date,
                gdf = gdf
            )

            self.centerlines.append(cl)

    def trim_centerlines_to_bounds(self):

        for cl in self.centerlines:
            cl.gdf = cl.trim_to_river_bounds(self.aoi)

    def merge_short_centerlines(self, kwargs):

        for cl in self.centerlines:
            cl.gdf = cl.merge_short_centerlines(**kwargs)

    def join_centerlines(self, kwargs):

        for cl in self.centerlines:
            cl.gdf = cl.join_cl_at_joints(**kwargs)

    def export_centerlines(self, file_type='geojson'):

        for cl in self.centerlines:
            path = os.path.join(self.outdir, f"centerlines_{cl.river_name}_{cl.date.strftime('%Y-%m-%d')}.{file_type}")
            cl.gdf.to_file(path)

    def add_swot(self, swot_path:str, date:datetime):

        swot = SWOT(
            river_name = self.name,
            date = date,
            path = swot_path,
            odir = self.outdir,
        )

        self.swot_obs.append(swot)

    def process_swot(self, dilate=5, engine='h5netcdf'):

        for i in range(len(self.masks)):

            # process the extraction mask for each mask
            self.masks[i].process_extraction_mask(dilate=dilate)

            # read in swot data and trim to extraction mask
            extraction_mask = self.masks[i].extraction_mask
            mask_transform = self.masks[i].transform
            mask_crs = self.masks[i].crs
            self.swot_obs[i].load_pixc(extraction_mask, mask_transform, mask_crs, engine=engine)


    def extract_water_levels(self):

        pass
        
