import os
import shutil

import geopandas as gpd
import numpy as np
import rasterio
import rioxarray
from pygeohydro import WBD
from pynhd import NHD
import py3dep
from skimage import morphology
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

def vectorize_stream(wbt, stream_binary_raster):
    """ stream raster should contain a single stream """
    thinned = morphology.thin(stream_binary_raster.data.astype(np.int8))
    thinned = thinned.astype(np.int8)

    with rasterio.open('temp.tif', 'w', driver='GTiff', height = thinned.shape[0],
                       width = thinned.shape[1], count=1, dtype=str(thinned.dtype),
                       crs=stream_binary_raster.rio.crs, transform=stream_binary_raster.rio.transform()) as dst:
        dst.write(thinned, 1)

    wbt.raster_to_vector_lines(os.path.abspath('temp.tif'), 'temp.shp')
    stream = gpd.read_file(os.path.join(wbt.work_dir, 'temp.shp'), crs=stream_binary_raster.rio.crs)

    os.remove('temp.tif')
    return stream

def get_stream_raster_from_hillslope(wbt, hillslope_raster):
    """ hillslope raster should contain a single hillslope """
    stream = hillslope_raster
    stream = stream.where(stream == 0)
    stream = stream.where(stream != 0, 1)
    return stream
