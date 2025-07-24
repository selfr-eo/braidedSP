![BraidedSP](images/braidedsp.png)

This repository holds the tools and code for the analysis presented within the pre-print [Franze et al. 2025](https://essopenarchive.org/doi/full/10.22541/essoar.174558891.16884958). The repo is designed to both provide the original notebooks used to conduct the analysis and create the figures used within the pre-print and also to organize the tools in a way that is easy to add your own data and process similar results for a river of your choosing.

## Publication archive
All notebooks, functions and figures related to the pre-print [Franze et al. 2025](https://essopenarchive.org/doi/full/10.22541/essoar.174558891.16884958) at the time of publishing can be found on the [publication folder]() within the publication branch. Find the README.md, within this folder within to learn more about the notebooks and tools.

> [!CAUTION]
> Tools and functionality on the main branch may have been updated, changed or improved since publication. If you are interested in the state of the code at the time of publishing make sure that you are viewing the code on the `publication branch`

## Why BraidedSP?
- The standard RiverSP product utilized simplified river centerlines that do not properly describe braided river networks, leading to errors and uncertainties in this product over braided rivers.
- This work utilized multispectral data from Sentinel-2 to define river masks and derive updated centerlines for braided rivers. The updated centerlines are then used to create a similar 'braidedSP' product that provides node aggregated elevations and widths along braided centerlines. See Figure 1 from [Preprint: Franze et al. 2025](https://essopenarchive.org/doi/full/10.22541/essoar.174558891.16884958).
- Understanding the differences water surface levels within various branches in braided rivers is crucial for proper hydraulic modeling and morphological studies. We found that there are significant differences between parallel channels in low flow conditions that are not seen by the RiverSP product. Furthermore comparison to high flows shows less difference between river slopes between such channels. See Figure 2 from [Preprint: Franze et al. 2025](https://essopenarchive.org/doi/full/10.22541/essoar.174558891.16884958).

![Figure 1](.\publication\figs\braidedSP_output.pdf)
![Figure 2](.\publication\figs\slopes_updated_07232025.pdf)


## BraidedSP Package and tools
A packaged version of the tools may be found within [braidedSP](). This packages organizes the tools and concepts implemented in a (hopefully) user friendly way to make processing your own braided river product easy.

### Getting started

#### Installation
> [!CAUTION]
> braidedSP has had limited testing but further developments may come. Installation via the following methods is not garunteed to work and may requre trouble shooting on your part.

While the package can be install directly from github using pip, there are a few caveats. The package relies upon GDAL which cannot be installed with a simple 'pip install gdal' command. We recommend installing gdal prior to package installation via your normal route or through conda-forge with the following command. We have not tested different gdal version but have pinned this specific gdal version to ensure compatability with the rasterio version pinned within the pyproject.toml

> Note: Ensure that `conda` or `mamba` are installed on your system before running the following command.
```sh
conda install gdal=3.10.2
```

To install the package directly from github:
> Note: Ensure that `pip` and `git` are installed on your system before running the following command.
```sh
pip install git+https://github.com/selfr-eo/braided_rivers
```

Alternatively, the pyproject.toml is built to support pixi. You can clone this directory and run the following command within the root directory of the folder.
> Note: Ensure that `pixi` is installed on your system before running the following command.
```sh
pixi install
```

#### Examples
It is recommended to check out our example workflow notebook ([full_workflow.ipynb]()) in the [examples]() folder which describes in detail what is needed to get started with the package. There are a few pieces of prior information needed, so no minimal working example is provided here.

#### Downloading SWOT PIXC data
Due to large size, example SWOT pixel cloud tiles are not included in this repository. To properly run the examples here, you should download the following tiles from wherever you get your SWOT data. We recommend using NASA's Earth Data search for quick downloads when you only need a few specific tiles: https://search.earthdata.nasa.gov/search. Choose project **SWOT**, collection **SWOT Level 2 Water Mask Pixel Cloud Data Product, Version C**, and paste the following for **Granule ID**.

Low flow tile (02-03-2024)
- SWOT_L2_HR_PIXC_010_258_112L_20240203T045957_20240203T050008_PIC0_01.nc

High flow tile (07-19-2024)
- SWOT_L2_HR_PIXC_018_258_112L_20240719T030033_20240719T030044_PIC0_01.nc

#### Processing and downloading S2 NDWI composites
Additionally, we do not provide the exact monthly Sentinel-2 NDWI composites used but provide a link to the snapshot code used to produce them. You can edit the geometry and dates to make a mask for your area or for the area in the examples:
https://code.earthengine.google.com/c2c5575a6ba1b2d728cdd5026689bb85





## Cite braidedSP

In case you use braidedSP in your research or work, it would be highly appreciated if you include a reference to our [Pre-Print](https://essopenarchive.org/doi/full/10.22541/essoar.174558891.16884958) in any kind of publication.

Cite as: Sarah Elizabeth Franze, Connor Chewning, Simon Jakob Köhn, Karina Nielsen. Braided rivers from SWOT: Water surface dynamics over multi-channel rivers. ESS Open Archive . April 25, 2025.
DOI: 10.22541/essoar.174558891.16884958/v1
<!-- 
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
``` -->