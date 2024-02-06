import os
import shutil

import geopandas as gpd
import numpy as np
import rioxarray

from valleys.utils import setup_wbt
from valleys.watershed import Watershed

# ------------ INPUTS ------------

dem_file =  "./sampledata/180600060101/dem.tif"
dem = rioxarray.open_rasterio(dem_file)
nhd_network_file = "./sampledata/180600060101/flowlines_mr.shp"
nhd_network = gpd.read_file(nhd_network_file)

wbt = setup_wbt(os.path.expanduser("~/opt/WBT/"), 
                os.path.abspath("./data/working/"))

# ------------ Watershed ------------
watershed = Watershed(dem, nhd_network, "./data/working/watershed/")
watershed.process_watershed(wbt)

# ------------ Valleys ------------
s_ids = np.unique(watershed.dataset.subbasins.values)
s_ids = s_ids[~np.isnan(s_ids)]
for s_id in s_ids:
    subbasin = watershed.clip_to_subbasin(s_id)

