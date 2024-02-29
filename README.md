# pyvalleys

This package is used to delineate the valley floor polygons of a river network.

# Setup

1. Install poetry
2. Install python 3.10+, recommend use pyenv and set project-specific python version
3. Install whiteboxtools and record the path to the directory containing the whitebox executable
4. run poetry install

For some reason, initializing the whiteboxtools class from 'whitebox' will 
download whiteboox tools. Even if you intend to directly pass the executable path to the instance.
This will only happen the first time and can be ignored.

# Usage

Input is a 10m DEM raster and a vector of the river network. 
Parameters for the algorithm are passed from a toml config file.
Output is a shapefile of the valley floor polygons as well as a directory with all the derived terrain attributes used for the algorithm.

This code can be used as a library or invoked from the command line.

Command line usage:
```
poetry run python -m pyvalleys <dem_path> <flowlines_path> <config_path> <wbt_exe_path> <terrain_attributes_dir> <ofile>
```

The workflow for extracting valleys can be seen in 'pyvalleys/workflow.py'.

# Config

An example config.toml file:

```
tolerance = 20 # for the simplification of the flowline
xs_spacing = 10 # spacing of cross-sections
xs_width = 100  # total width is twice this value
xs_point_spacing = 10 # spacing of points along cross-sections
quantile = 0.75 # sets the hand threshold
buffer = 2 # add to hand threshold 
slope_threshold = 14.1 # a little more than 25% slope

peak_threshold = 0.002
bp_slope_threshold = 20
```
