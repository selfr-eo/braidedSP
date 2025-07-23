![BraidedSP](images/braidedsp.png)

This repository holds the tools and code for the analysis presented within the paper [Franze et al. 2025](link_to_DOI). The repo is designed to both provide the original notebooks used to conduct the analysis and create the figures used within the paper and also to organize the tools in a way that is easy to add your own data and process similar results for a river of your choosing.

## Publication archive
All notebooks, functions and figures related to the paper [Franze et al. 2025](link_to_DOI) at the time of publishing can be found on the 'publication' branch. Find the README.md, within the 'publication folder' to learn more about the notebooks and tools.

NOTE: Tools and functionality on the main branch may have been updated, changed or improved since publication. If you are interested in the state of the code at the time of publishing make sure that you are viewing the code on the 'publication branch'

## BraidedSP Package and tools
A packaged version of the tools may be found within [braidedSP](). This packages organizes the tools and concepts implemented in a (hopefully) user friendly way to make processing your own braided river product easy.

### Getting started

#### Installation
> Note: Ensure that `pip` and `git` are installed on your system before running the following command.
```sh
pip install git+https://github.com/selfr-eo/braided_rivers
```

#### Examples
It is recommended to check out our example workflow notebook ([full_workflow.ipynb]()) in the [examples]() folder which describes in detail what is needed to get started with the package. There are a few pieces of prior information needed, so no minimal working example is provided here.



##### Downloading SWOT PIXC data
Due to large size, example SWOT pixel cloud tiles are not included in this repository. To properly run the examples here, you should download the following tiles from wherever you get your SWOT data. We recommend using NASA's Earth Data search for quick downloads when you only need a few specific tiles: https://search.earthdata.nasa.gov/search. Choose project **SWOT**, collection **SWOT Level 2 Water Mask Pixel Cloud Data Product, Version C**, and paste the following for **Granule ID**.

Low flow tile (02-03-2024)
- SWOT_L2_HR_PIXC_010_258_112L_20240203T045957_20240203T050008_PIC0_01.nc

High flow tile (07-19-2024)
- SWOT_L2_HR_PIXC_018_258_112L_20240719T030033_20240719T030044_PIC0_01.nc

##### Processing and downloading S2 NDWI composites
Additionally, we do not provide the exact monthly Sentinel-2 NDWI composites used but provide a link to the snapshot code used to produce them. You can edit the geometry and dates to make a mask for your area or for the area in the examples:
https://code.earthengine.google.com/c2c5575a6ba1b2d728cdd5026689bb85





## Cite braidedSP

In case you use braidedSP in your research or work, it would be highly appreciated if you include a reference to our [Paper](link_to_DOI) in any kind of publication.

```bibtex
@article{franze2025,
  title = {},
  author = {},
  journal = {},
  publisher = {},
  year = {},
  volume = {},
  number = {},
  pages = {},
  doi = {},
  url = {},
}
```