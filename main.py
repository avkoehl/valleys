import os
import shutil
import glob

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

def delineate_valleys(watershed, odir):
    if os.path.isdir(odir):
        shutil.rmtree(odir)
    os.makedirs(odir)

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
    bps = pd.concat(bps)

    valley_floors_df.to_file(os.path.join(odir, "valleys.shp"))
    bps.to_file(os.path.join(odir, "breakpoints.shp"))

def full_workflow(dem_file, nhd_network_file, wbt, odir):
    dem = rioxarray.open_rasterio(dem_file)
    nhd_network = gpd.read_file(nhd_network_file)

    watershed = Watershed(dem, nhd_network, os.path.join(odir))
    watershed.process_watershed(wbt)
    delineate_valleys(watershed, os.path.join(odir, 'results/'))

wbt = setup_wbt('~/opt/WBT/', './data/working/')

huc12s = ['180400060301', '180600060101', '180701060201']
for huc12 in huc12s:
    dem_file = f"./sampledata/{huc12}/dem.tif"
    nhd_network_file = f"./sampledata/{huc12}/flowlines_mr.shp"
    odir = f"./data/{huc12}/"
    full_workflow(dem_file, nhd_network_file, wbt, odir)
