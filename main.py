import os
import shutil

import geopandas as gpd
import rioxarray
from shapely.geometry import box
import pandas as pd

from valleys.utils import setup_wbt
from valleys.cross_section import get_cross_section_points
from valleys.terrain import flow_accumulation_workflow
from valleys.terrain import align_to_terrain
from valleys.terrain import segment_subbasins
from valleys.terrain import clip_to_subbasin
from valleys.valley import add_attributes_to_xs
from valleys.valley import get_break_points

# ------------ INPUTS ------------

dem_file =  "./sampledata/180600060101/dem.tif"
nhd_network_file = "./sampledata/180600060101/flowlines_mr.shp"

wbt = setup_wbt(os.path.expanduser("~/opt/WBT/"), 
                os.path.abspath("./data/working/"))

# ------------ SEGMENT ------------

nhd_network = gpd.read_file(nhd_network_file)

flow_acc_files = flow_accumulation_workflow(wbt, dem_file)
stream_files = align_to_terrain(wbt, nhd_network, flow_acc_files['flow_accum'], flow_acc_files['flow_dir'])
subbasin_files = segment_subbasins(wbt, stream_files['streams'], flow_acc_files['flow_dir'])

# load rasters
subbasins_f = rioxarray.open_rasterio(subbasin_files['subbasins'], masked=True).squeeze()
flowlines_f = gpd.read_file(stream_files['streams_shp'])
dem_f = rioxarray.open_rasterio(flow_acc_files['conditioned_dem'], masked=True).squeeze()
hillslopes_f = rioxarray.open_rasterio(subbasin_files['hillslopes'], masked=True).squeeze()
flow_dir_f = rioxarray.open_rasterio(flow_acc_files['flow_dir'], masked=True).squeeze()
stream_f = rioxarray.open_rasterio(stream_files['streams'], masked=True).squeeze()

# ------------ VALLEY DELINEATION ------------
count = 0
for subbasin_id in pd.Series(subbasins.values.flatten()).dropna().unique():
    print(subbasin_id)
    flowline = flowlines.loc[flowlines['STRM_VAL'] == subbasin_id]['geometry'].iloc[0]
    dem, hillslope, flow_dir, stream = clip_to_subbasin(dem_f, hillslopes_f, flow_dir_f, stream_f, subbasins_f, subbasin_id)

    points = get_cross_section_points(flowline, xs_spacing=50, xs_width=500, xs_point_spacing=10, tolerance=20)
    points = add_attributes_to_xs(wbt, points, dem, stream, flow_dir, hillslope)

    break_points_df = get_break_points(points)
    if count == 5:
        break
    count = count + 1
