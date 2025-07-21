![BraidedSP](images/braidedsp.png)

This packages holds the tools and code for reproducing the figures and analysis presented within the paper [Franze et al. 2025](link_to_DOI). It has been structured in a way that it should be easy to add your own data and process similar results for a river of your choosing.

## Cite braidedSP

In case you use baidedSP in your research or work, it would be highly appreciated if you include a reference to our [Paper](link_to_DOI) in any kind of publication.

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

## Installation
TODO

## Outline of raw notebooks:

- gen_centerlines.ipynb (and gen_centerlines_niger.ipynb setup specifically for Niger) - notebook for generating centerlines from S2 images. Adjust paths to the appropriate data.

- gen_BraidedSP.ipnb - notebook for generating the braidedSP product over the created centerlines.

- view_BraidedSP - notebook for viewing the processed SWOT data with respect to the generated centerlines. Paper figures plotted here.

## Downloading SWOT PIXC data
Due to large size, example SWOT pixel cloud tiles are not included in this repository. To properly run the examples here, you should download the following tiles from wherever you get your SWOT data. We recommend using NASA's Earth Data search for quick downloads when you only need a few specific tiles: https://search.earthdata.nasa.gov/search. Choose project **SWOT**, collection **SWOT Level 2 Water Mask Pixel Cloud Data Product, Version C**, and paste the following for **Granule ID**.

Low flow tile (02-03-2024)
- SWOT_L2_HR_PIXC_010_258_112L_20240203T045957_20240203T050008_PIC0_01.nc

High flow tile (07-19-2024)
- SWOT_L2_HR_PIXC_018_258_112L_20240719T030033_20240719T030044_PIC0_01.nc