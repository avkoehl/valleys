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

    bps = []
    results_list = []
    for sid in watershed.get_subbasin_ids():
        results = {}
        print(sid)
        subbasin_data, flowline = watershed.clip_to_subbasin(sid)
        subbasin_data = prep_dataset(subbasin_data)
    
        subbasin = Subbasin(subbasin_data, flowline, sid)
        subbasin.sample_cross_section_points()
        subbasin.find_breakpoints()

        results['subbasin_id'] = sid

        # options:
        #  quantile = 0.5
        #  quantile = 0.5 slope_threshold = 25
        #  quantile = 0.5 slope_threshold = 25 buffer = 5
        subbasin.delineate_valley_floor(quantile=0.5, slope_threshold=None, buffer=0)
        results['q50_s0_b0'] = subbasin.valley_floor_polygon

        subbasin.delineate_valley_floor(quantile=0.5, slope_threshold=25, buffer=0)
        results['q50_s25_b0'] = subbasin.valley_floor_polygon

        subbasin.delineate_valley_floor(quantile=0.5, slope_threshold=25, buffer=5)
        results['q50_s25_b5'] = subbasin.valley_floor_polygon


        #  quantile = 0.75  
        #  quantile = 0.75 slope_threshold = 25
        #  quantile = 0.75 slope_threshold = 25 buffer = 5
        subbasin.delineate_valley_floor(quantile=0.75, slope_threshold=None, buffer=0)
        results['q75_s0_b0'] = subbasin.valley_floor_polygon

        subbasin.delineate_valley_floor(quantile=0.75, slope_threshold=25, buffer=0)
        results['q75_s25_b0'] = subbasin.valley_floor_polygon

        subbasin.delineate_valley_floor(quantile=0.75, slope_threshold=25, buffer=5)
        results['q75_s25_b5'] = subbasin.valley_floor_polygon

        #  quantile = 0.95  
        #  quantile = 0.95  slope_threshold = 25
        #  quantile = 0.95  slope_threshold = 25 buffer = 5
        subbasin.delineate_valley_floor(quantile=0.95, slope_threshold=None, buffer=0)
        results['q95_s0_b0'] = subbasin.valley_floor_polygon

        subbasin.delineate_valley_floor(quantile=0.95, slope_threshold=25, buffer=0)
        results['q95_s25_b0'] = subbasin.valley_floor_polygon

        subbasin.delineate_valley_floor(quantile=0.95, slope_threshold=25, buffer=5)
        results['q95_s25_b5'] = subbasin.valley_floor_polygon

        results_list.append(results)
        bps.append(subbasin.break_points_df)

    df = pd.DataFrame.from_records(results_list)
    q50_s0_b0 = gpd.GeoDataFrame(df[['subbasin_id', 'q50_s0_b0']], geometry='q50_s0_b0',  crs=3310)
    q50_s25_b0 = gpd.GeoDataFrame(df[['subbasin_id', 'q50_s25_b0']], geometry='q50_s25_b0',  crs=3310)
    q50_s25_b5 = gpd.GeoDataFrame(df[['subbasin_id', 'q50_s25_b5']], geometry='q50_s25_b5',  crs=3310)
    q75_s0_b0 = gpd.GeoDataFrame(df[['subbasin_id', 'q75_s0_b0']], geometry='q75_s0_b0',  crs=3310)
    q75_s25_b0 = gpd.GeoDataFrame(df[['subbasin_id', 'q75_s25_b0']], geometry='q75_s25_b0',  crs=3310)
    q75_s25_b5 = gpd.GeoDataFrame(df[['subbasin_id', 'q75_s25_b5']], geometry='q75_s25_b5',  crs=3310)
    q95_s0_b0 = gpd.GeoDataFrame(df[['subbasin_id', 'q95_s0_b0']], geometry='q95_s0_b0',  crs=3310)
    q95_s25_b0 = gpd.GeoDataFrame(df[['subbasin_id', 'q95_s25_b0']], geometry='q95_s25_b0',  crs=3310)
    q95_s25_b5 = gpd.GeoDataFrame(df[['subbasin_id', 'q95_s25_b5']], geometry='q95_s25_b5',  crs=3310)
    bps = pd.concat(bps)

    q50_s0_b0.to_file(os.path.join(odir, "q50_s0_b0.shp"))
    q50_s25_b0.to_file(os.path.join(odir, "q50_s25_b0.shp"))
    q50_s25_b5.to_file(os.path.join(odir, "q50_s25_b5.shp"))
    q75_s0_b0.to_file(os.path.join(odir, "q75_s0_b0.shp"))
    q75_s25_b0.to_file(os.path.join(odir, "q75_s25_b0.shp"))
    q75_s25_b5.to_file(os.path.join(odir, "q75_s25_b5.shp"))
    q95_s0_b0.to_file(os.path.join(odir, "q95_s0_b0.shp"))
    q95_s25_b0.to_file(os.path.join(odir, "q95_s25_b0.shp"))
    q95_s25_b5.to_file(os.path.join(odir, "q95_s25_b5.shp"))

    bps.to_file(os.path.join(odir, "breakpoints.shp"))
    watershed.dataset['dem'].rio.to_raster(os.path.join(odir, "dem.tif"))
    watershed.dataset['flowpaths_identified'].rio.to_raster(os.path.join(odir, "streams.tif"))

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
