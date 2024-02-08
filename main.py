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

    valley_floors_1 = []
    valley_floors_2 = []
    valley_floors_3 = []
    bps = []
    for sid in watershed.get_subbasin_ids():
        print(sid)
        subbasin_data, flowline = watershed.clip_to_subbasin(sid)
        subbasin_data = prep_dataset(subbasin_data)
    
        subbasin = Subbasin(subbasin_data, flowline, sid)
        subbasin.sample_cross_section_points()
        subbasin.find_breakpoints()
        subbasin.determine_hand_threshold()

        # option 1: use hand threshold
        subbasin.delineate_valley_floor(slope_threshold=None, hand_buffer=0)
        valley_floors_1.append((sid, subbasin.valley_floor_polygon))
        # option 2: use hand threshold and slope threshold
        subbasin.delineate_valley_floor(slope_threshold=25, hand_buffer=0)
        valley_floors_2.append((sid, subbasin.valley_floor_polygon))
        # option 3: use hand threshold with buffer and slope threshold
        subbasin.delineate_valley_floor(slope_threshold=25, hand_buffer=10)
        valley_floors_3.append((sid, subbasin.valley_floor_polygon))
    
        bps.append(subbasin.break_points_df)
    
    valley_floors_df_1 = gpd.GeoDataFrame(valley_floors_1, columns=['subbasin_id', 'geometry'], crs=3310)
    valley_floors_df_2 = gpd.GeoDataFrame(valley_floors_2, columns=['subbasin_id', 'geometry'], crs=3310)
    valley_floors_df_3 = gpd.GeoDataFrame(valley_floors_3, columns=['subbasin_id', 'geometry'], crs=3310)
    bps = pd.concat(bps)

    valley_floors_df_1.to_file(os.path.join(odir, "valleys_hand.shp"))
    valley_floors_df_2.to_file(os.path.join(odir, "valleys_slope.shp"))
    valley_floors_df_3.to_file(os.path.join(odir, "valleys_buffer_and_slope.shp"))
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
