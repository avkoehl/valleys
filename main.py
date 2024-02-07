import os
import shutil

import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray

from valleys.utils import setup_wbt
from valleys.watershed import Watershed
from valleys.subbasin import Subbasin

def prep_dataset(dataset):
    mapping = {
        'conditioned_dem': 'elevation',
        'flowpaths_identified': 'strm_val'
    }
    # rename bands to match mapping
    dataset = dataset.rename(mapping)
    keys = ['elevation', 'slope', 'curvature', 'strm_val', 'hillslopes', 'flow_dir', 'hand']
    dataset = dataset[keys]
    return dataset

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
valley_floors = []
bps = []
for sid in watershed.get_subbasin_ids():
    print(sid)
    subbasin_data, flowline = watershed.clip_to_subbasin(sid)
    subbasin_data = prep_dataset(subbasin_data)

    subbasin = Subbasin(subbasin_data, flowline, sid)
    subbasin.valley_floor_by_breakpoints_full_workflow()

    valley_floors.append((sid, subbasin.valley_floor_polygon))
    bps.append(subbasin.break_points_df)

valley_floors_df = gpd.GeoDataFrame(valley_floors, columns=['subbasin_id', 'geometry'], crs=3310)
valley_floors_df.to_file("./data/outputs/valleys.shp")

pd.concat(bps).to_file("breakpoints.shp")

