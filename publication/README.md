
# Publication archive for [Franze et al. 2025](link_to_DOI)

This directory holds all of the notebooks and scripts for the analysis conducted within [Franze et al. 2025](link_to_DOI). While this folder should not be edited in future editions of the repository, to ensure that you are viewing the tools related to the paper at the time of publishing please make sure that you are viewing this on the 'publication' branch.

Note: the scripts provided here are to aid in reproduceability of results, but are not organized or written in a way that is necessarily easy to run. Thus we have provided the tools in a packaged format. If you are more interested in running the workflow on your own area of interest, information on the package can be found on within the main directory of this repo. [braidedSP]().

Additional note: Due to the large file size of the SWOT PIXC data and Sentinel-2 composites, it is unefeasible to share the exact files ingested within these scripts. We have instead provided instructions on how to extract the same SWOT PIXC data as we used directly from NASA's Earth Data. See [Downloading SWOT PIXC data](####Downloading-SWOT-PIXC-data).


## Directory structure

- [notebooks](): Holds the notebooks used for analysis of data.
    - [gen_centerlines.ipynb]() - Notebook used to generate braided centerlines from monthly sentinel-2 NDWI composites
    - [gen_BraidedSP.ipynb]() - Notebook used to extract SWOT PIXC data to braidedSP product
    - [view_BraidedSP.ipynb]() - Notebook used to compare braidedSP product with SWOT RiverSP product
    - [width_analysis.ipynb]() - Additional notebook used to analyse widths extracted BraidedSP product
- [tools](): Additional tools used within the notebooks
    - [procBraided.py]()
    - [skeletonize_func.py]()
-[figs](): Some of the raw figures produced by the scripts that are presented within the paper.


## Downloading SWOT PIXC data
Due to large size, example SWOT pixel cloud tiles are not included in this repository. To properly run the examples here, you should download the following tiles from wherever you get your SWOT data. We recommend using NASA's Earth Data search for quick downloads when you only need a few specific tiles: https://search.earthdata.nasa.gov/search. Choose project **SWOT**, collection **SWOT Level 2 Water Mask Pixel Cloud Data Product, Version C**, and paste the following for **Granule ID**.

Low flow tile (02-03-2024)
- SWOT_L2_HR_PIXC_010_258_112L_20240203T045957_20240203T050008_PIC0_01.nc

High flow tile (07-19-2024)
- SWOT_L2_HR_PIXC_018_258_112L_20240719T030033_20240719T030044_PIC0_01.nc

## Processing and downloading S2 NDWI composites
Additionally, we do not provide the exact monthly Sentinel-2 NDWI composites used but provide a link to the snapshot code used to produce them. You can edit the geometry and dates to make a mask for your area or for the area in the examples:
https://code.earthengine.google.com/c2c5575a6ba1b2d728cdd5026689bb85

