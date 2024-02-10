import os
import shutil
import glob

import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray

from valleys.utils import setup_wbt
from valleys.workflow import prep_dataset
from valleys.workflow import delinate_valleys
from valleys.watershed import Watershed
from valleys.subbasin import Subbasin

dem_file = './data/1803001003/dem.tif'
nhd_network_file = './data/1803001003/flowlines_mr.shp'
wbt = setup_wbt('~/opt/WBT/', './data/working/')
odir = './data/valleys_1803001003/'

if os.path.exists(odir):
    shutil.rmtree(odir)
os.makedirs(odir)

dem = rioxarray.open_rasterio(dem_file)
nhd_network = gpd.read_file(nhd_network_file)

watershed = Watershed(dem, nhd_network, wbt.work_dir)
watershed.process_watershed(wbt)
valleys = delinate_valleys(watershed)

valleys.to_file(os.path.join(odir, 'valley_floor.shp'))
watershed.dataset.dem.rio.to_raster(os.path.join(odir, 'dem.tif'))
watershed.dataset.flowpaths_identified.rio.to_raster(os.path.join(odir, 'flowpaths.tif'))
