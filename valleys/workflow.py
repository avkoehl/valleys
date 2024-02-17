import datetime
import os
import shutil

import geopandas as gpd
import rioxarray
import toml

from valleys.subbasin import Subbasin
from valleys.utils import setup_wbt
from valleys.watershed import Watershed

def full_workflow(config):
    # store config files as toml
    # required keys in config
    required_keys = ['wbt_path', 'dem_path', 'nhd_mr_path', 'output_dir', 'params']
    required_keys_params = ['tolerance', 'xs_spacing', 'xs_width', 'xs_point_spacing', 'quantile', 'buffer', 'slope_threshold']
    for key in required_keys:
        if key not in config:
            raise ValueError(f'config is missing required key: {key}')
    for key in required_keys_params:
        if key not in config['params']:
            raise ValueError(f'config is missing required key: {key}')

    WBT_WORK_DIR = os.path.join(config['output_dir'], 'wbt_work_dir')

    # setup directories
    output_dir = config['output_dir']
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    os.makedirs(WBT_WORK_DIR)

    # setup wbt
    wbt = setup_wbt(config['wbt_path'], WBT_WORK_DIR)

    # load data
    dem = rioxarray.open_rasterio(config['dem_path'])
    nhd_mr = gpd.read_file(config['nhd_mr_path'])

    # create watershed and extract valleys
    watershed = Watershed(dem, nhd_mr, WBT_WORK_DIR)
    watershed.process_watershed(wbt)
    valleys = delineate_valleys(watershed, **config['params'])

    # save valleys
    valleys.to_file(os.path.join(output_dir, 'valley_floors.shp'))
    watershed.dataset.dem.rio.to_raster(os.path.join(output_dir, 'dem.tif'))
    watershed.dataset.flowpaths_identified.rio.to_raster(os.path.join(output_dir, 'flowpaths.tif'))

    # save config as toml to output dir
    # get current date
    date = datetime.datetime.now().strftime("%Y-%m-%d")
    with open(os.path.join(output_dir, f'{data}_config.toml'), 'w') as f:
        toml.dump(config, f)
    return 

def delineate_valleys(watershed,
                     tolerance=20,
                     xs_spacing=50,
                     xs_width = 500,
                     xs_point_spacing = 10,
                     quantile = 0.75,
                     buffer = 3,
                     slope_threshold = 20):

    results = []
    redo = []
    for sid in watershed.get_subbasin_ids():
        subbasin_data, flowline = watershed.clip_to_subbasin(sid)
        subbasin_data = prep_dataset(subbasin_data)

        subbasin = Subbasin(subbasin_data, flowline, sid)
        subbasin.valley_floor_by_breakpoints_full_workflow(
                tolerance=tolerance,
                xs_spacing=xs_spacing,
                xs_width = xs_width,
                xs_point_spacing = xs_point_spacing,
                quantile = quantile,
                buffer = buffer,
                slope_threshold = slope_threshold)
        poly = subbasin.valley_floor_polygon
        threshold = subbasin.hand_threshold
        if threshold is None:
            redo.append(sid)

        if threshold is not None:
            results.append((sid, poly, threshold, quantile, buffer, slope_threshold))
    hand_thresholds = [r[2] for r in results]
    mean = sum(hand_thresholds) / len(hand_thresholds)
    for sid in redo:
        subbasin_data, flowline = watershed.clip_to_subbasin(sid)
        subbasin_data = prep_dataset(subbasin_data)
        subbasin = Subbasin(subbasin_data, flowline, sid)
        subbasin.delineate_valley_floor(buffer=buffer, slope_threshold=slope_threshold, overwrite_hand=mean)
        results.append((sid, poly, threshold, quantile, buffer, slope_threshold))

    df = gpd.GeoDataFrame(results, columns=['ID', 'floor', 'HAND', 'quantile', 'buffer', 'max_slope'], geometry='floor')
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
