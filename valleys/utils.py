import os
import shutil

import geopandas as gpd
import rioxarray
import whitebox

def setup_wbt(whitebox_dir, working_dir):
    wbt = whitebox.WhiteboxTools()
    wbt.set_whitebox_dir(os.path.expanduser(whitebox_dir))
    
    working_directory = os.path.abspath(working_dir)
    if os.path.exists(working_directory):
            shutil.rmtree(working_directory)
    os.mkdir(working_directory)
    wbt.set_working_dir(os.path.abspath(working_directory))
    wbt.set_verbose_mode(False)
    return wbt

def chomp_raster(raster):
    raster = raster.dropna(dim='x', how='all')
    raster = raster.dropna(dim='y', how='all')
    return raster

def get_extent(raster):
    """
     raster: xarray.DataArray
     raster.bounds() is in format (left, bottom, right, top)
     returns: tuple of (left, bottom, right, top)
    """
    bounds = raster.rio.bounds()
    return (bounds[0], bounds[2], bounds[1], bounds[3])

