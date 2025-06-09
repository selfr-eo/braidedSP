# Typing imports
from dataclasses import dataclass
from datetime import datetime

# general imports
from tqdm import tqdm
import geopandas as gpd

# import project specific objects
from braidedSP.mask import Mask
from braidedSP.centerline import Centerline


@dataclass
class River:

    name: str
    aoi: gpd.GeoDataFrame
    outdir: str

    def __post_init__(self):

        self.masks = []
        self.centerlines = []

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

        # cycle through and process masks
        for mask in tqdm(self.masks, desc='Processing masks', leave=True, disable=not show_progress):
            mask.process(**kwargs)

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