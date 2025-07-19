![BraidedSP]('./images/braidedsp.png')

# Braided Rivers SP

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

