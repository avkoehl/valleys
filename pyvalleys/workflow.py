import datetime
import os
import shutil

import geopandas as gpd
import rioxarray
import toml
import whitebox

from pyvalleys.subbasin import Subbasin
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

def valley_floors(dem_file, flowlines_file, config_file, wbt_path, terrain_dir, ofile):
    dem, flowlines, config, wbt, terrain_dir, ofile = setup(dem_file, flowlines_file, config_file, wbt_path, terrain_dir, ofile)
    watershed = Watershed(dem, flowlines, terrain_dir)
    watershed.process_watershed(dem)
    valleys = delineate_valleys(watershed, **config)
    valleys['date'] = datetime.datetime.now().strftime("%Y-%m-%d")
    valleys['config'] = config
    valleys['version'] = valleys.__version__
    valleys['wbt_version'] = wbt.version()
    valleys.to_file(ofile)

def setup(dem_file, flowlines_file, config_file, wbt_path, terrain_dir, ofile):
    # convert to absolute paths
    dem_file = os.path.abspath(dem_file)
    flowlines_file = os.path.abspath(flowlines_file)
    config_file = os.path.abspath(config_file)
    wbt_path = os.path.abspath(wbt_path)
    terrain_dir = os.path.abspath(terrain_dir)
    ofile = os.path.abspath(ofile)

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
                     tolerance=20,
                     xs_spacing=50,
                     xs_width = 500,
                     xs_point_spacing = 10,
                     quantile = 0.75,
                     buffer = 3,
                     slope_threshold = 20,
                     peak_threshold = 0.002,
                     bp_slope_threshold = 20):


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
                slope_threshold = slope_threshold,
                peak_threshold = peak_threshold,
                bp_slope_threshold = bp_slope_threshold
                )

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

    df = gpd.GeoDataFrame(results, columns=['ID', 'floor', 'HAND', 'quantile', 'buffer', 'max_slope'], geometry='floor', crs=watershed.dataset.dem.rio.crs)
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
