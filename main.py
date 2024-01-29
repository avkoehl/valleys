import os
import shutil

import geopandas as gpd

from valleys.fun import align_to_terrain, compute_terrain_rasters, network_cross_section_points
from valleys.utils import setup_wbt

dem_file =  "./sampledata/180600060101/dem.tif"
nhd_network_file = "./sampledata/180600060101/flowlines_mr.shp"
nhd_network = gpd.read_file(nhd_network_file)

wbt_exe_path = os.path.expanduser("~/opt/WBT/")
working_dir = os.path.abspath("./data/working/")

if os.path.isdir(working_dir):
    shutil.rmtree(working_dir)
os.makedirs(working_dir)

wbt = setup_wbt(wbt_exe_path, working_dir)

# copy dem and nhd to working dir so that all the outputs are in the same place

paths = align_to_terrain(wbt, nhd_network, dem_file)
flowlines = gpd.read_file(paths['flow_paths_shp'])

rasters = compute_terrain_rasters(wbt, paths['conditioned_dem'],
                                  paths['flow_path'], paths['flow_dir'])

points = network_cross_section_points(flowlines)

for stream_id in points.stream_id.unique():
    stream_points = points.loc[points.stream_id == stream_id]
    #valley = delineate_valley(stream_id, stream_points, rasters)
