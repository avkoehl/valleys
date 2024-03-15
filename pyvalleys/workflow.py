import datetime
import os
import shutil
import sys
import json

import geopandas as gpd
import pandas as pd
import rioxarray
import toml
import whitebox

from pyvalleys.valley_extractor import ValleyExtractor
from pyvalleys.watershed import Watershed

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

def valley_floors(dem_file, flowlines_file, config_file, wbt_path, terrain_dir, ofile, debug_file=None):
    dem, flowlines, config, wbt, terrain_dir, ofile = setup(dem_file, flowlines_file, config_file, wbt_path, terrain_dir, ofile)
    watershed = Watershed(dem, flowlines, terrain_dir)
    watershed.process_watershed(wbt)

    if debug_file is not None:
        valleys, breakpoints = delineate_valleys(watershed, **config, debug_flag=True)
    else:
        valleys = delineate_valleys(watershed, **config)
    valleys['date'] = datetime.datetime.now().strftime("%Y-%m-%d")
    valleys['config'] = json.dumps(config)
    valleys.to_file(ofile)

    if debug_file is not None:
        breakpoints.to_file(debug_file)

def convert_to_absolute_path(filename):
    return os.path.abspath(os.path.expanduser(filename))

def setup(dem_file, flowlines_file, config_file, wbt_path, terrain_dir, ofile):
    # convert to absolute paths
    dem_file = convert_to_absolute_path(dem_file)
    flowlines_file = convert_to_absolute_path(flowlines_file)
    config_file = convert_to_absolute_path(config_file)
    wbt_path = convert_to_absolute_path(wbt_path)
    terrain_dir = convert_to_absolute_path(terrain_dir)
    ofile = convert_to_absolute_path(ofile)

    # make sure dem_file and flowlines_file exist
    if not os.path.exists(dem_file):
        sys.exit(f"DEM file {dem_file} does not exist")
    if not os.path.exists(flowlines_file):
        sys.exit(f"Flowlines file {flowlines_file} does not exist")

    # setup directories and whiteboxtools
    if not os.path.exists(terrain_dir):
        os.makedirs(terrain_dir)

    # get directory from ofile
    ofile_dir = os.path.dirname(ofile)
    if not os.path.exists(ofile_dir):
        os.makedirs(ofile_dir)

    wbt = setup_wbt(wbt_path, terrain_dir)

    dem = rioxarray.open_rasterio(dem_file)
    flowlines = gpd.read_file(flowlines_file)

    # parse config_file
    config = toml.load(config_file)
    required_keys = ['tolerance', 'xs_spacing', 'xs_width', 'xs_point_spacing', 'quantile', 'buffer', 'slope_threshold', 'peak_threshold', 'bp_slope_threshold']
    for key in required_keys:
        if key not in config:
            raise ValueError(f'config is missing required key: {key}')

    return dem, flowlines, config, wbt, terrain_dir, ofile


def delineate_valleys(watershed,
                      debug_flag=False,
                      smooth=True,
                     tolerance=20,
                     xs_spacing=50,
                     xs_width = 500,
                     xs_point_spacing = 10,
                     quantile = 0.75,
                     buffer = 3,
                     slope_threshold = 20,
                     peak_threshold = 0.002,
                     bp_slope_threshold = 20):


    if debug_flag:
        debug = [] # list of breakpoint dfs
    results = []
    for sid in watershed.get_subbasin_ids():
        print(sid)
        subbasin_data, flowline = watershed.clip_to_subbasin(sid)
        subbasin_data = prep_dataset(subbasin_data)

        ve = ValleyExtractor(subbasin_data, flowline, sid)
        ve.run(
                tolerance=tolerance,
                smooth= smooth,
                xs_spacing=xs_spacing,
                xs_width = xs_width,
                xs_point_spacing = xs_point_spacing,
                quantile = quantile,
                buffer = buffer,
                slope_threshold = slope_threshold,
                peak_threshold = peak_threshold,
                bp_slope_threshold = bp_slope_threshold
                )

        poly = ve.valley_floor_polygon
        threshold = ve.hand_threshold
        results.append((sid, poly, threshold, quantile, buffer, slope_threshold))

        if debug_flag:
            bps = ve.break_points_df
            bps['Subbasin_ID'] = sid
            debug.append(bps)

    # HOW TO HANDLE CASES WHERE NO HAND THRESHOLD?
    df = gpd.GeoDataFrame(results, columns=['ID', 'floor', 'HAND', 'quantile', 'buffer', 'max_slope'], geometry='floor', crs=watershed.dataset.dem.rio.crs)

    if debug_flag:
        bps = pd.concat(debug)
        return df, bps

    return df

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
