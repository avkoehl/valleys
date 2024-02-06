import os
import shutil

import geopandas as gpd
import numpy as np
import rioxarray

from valleys.utils import setup_wbt
from valleys.watershed import Watershed
from valleys.subbasin import Subbasin

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
    subbasin_data, flowline = watershed.clip_to_subbasin(s_id)
    mapping = {
            'conditioned_dem': 'elevation',
            'flowpaths_identified': 'strm_val'
            }
    # rename bands to match mapping
    subbasin_data = subbasin_data.rename(mapping)
    keys = ['elevation', 'slope', 'curvature', 'strm_val', 'hillslopes', 'flow_dir', 'hand']
    subbasin_data = subbasin_data[keys]
    subbasin = Subbasin(subbasin_data, flowline, s_id)

    subbasin.sample_cross_section_points()
    subbasin.find_breakpoints()
    subbasin.determine_hand_threshold()
    subbasin.delineate_valley_floor()
