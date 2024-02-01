import os
import shutil

import geopandas as gpd
import rioxarray
from shapely.geometry import box
import pandas as pd

from valleys.utils import setup_wbt
from valleys.utils import chomp_raster
from valleys.cross_section import get_cross_section_points
from valleys.fun import flow_accumulation_workflow
from valleys.fun import align_to_terrain
from valleys.fun import segment_subbasins

# ------------ INPUTS ------------

dem_file =  "./sampledata/180600060101/dem.tif"
nhd_network_file = "./sampledata/180600060101/flowlines_mr.shp"

wbt_exe_path = os.path.expanduser("~/opt/WBT/")
working_dir = os.path.abspath("./data/working/")

if os.path.isdir(working_dir):
    shutil.rmtree(working_dir)
os.makedirs(working_dir)

wbt = setup_wbt(wbt_exe_path, working_dir)

# ------------ WORKFLOW ------------

nhd_network = gpd.read_file(nhd_network_file)

# flow accumulation workflow 
flow_acc_files = flow_accumulation_workflow(wbt, dem_file)

# align nhd to terrain
stream_files = align_to_terrain(wbt, nhd_network, flow_acc_files['flow_accum'], flow_acc_files['flow_dir'])

# segment basin
subbasin_files = segment_subbasins(wbt, stream_files['streams'], flow_acc_files['flow_dir'])

# load subbasins raster and flowlines
subbasins = rioxarray.open_rasterio(subbasin_files['subbasins'], masked=True).squeeze()
flowlines = gpd.read_file(stream_files['streams_shp'])

dem = rioxarray.open_rasterio(flow_acc_files['conditioned_dem'], masked=True).squeeze()
hillslopes = rioxarray.open_rasterio(subbasin_files['hillslopes'], masked=True).squeeze()
flow_dir = rioxarray.open_rasterio(flow_acc_files['flow_dir'], masked=True).squeeze()
stream = rioxarray.open_rasterio(stream_files['streams'], masked=True).squeeze()

for row in flowlines.iterrows():
    subbasin_id = row['STRM_VAL']
    flowline = row['geometry']

    # clip dem, stream, and flow dir and hillslope to subbasin
    dem_clipped = dem.where(subbasins == subbasin_id)
    dem_clipped = chomp_raster(dem_clipped)
    
    hillslope_clipped = hillslopes.where(subbasins == subbasin_id)
    hillslope_clipped = chomp_raster(hillslope_clipped)

    flow_dir_clipped = flow_dir.where(subbasins == subbasin_id)
    flow_dir_clipped = chomp_raster(flow_dir_clipped)

    stream_clipped = stream.where(stream == subbasin_id)
    stream_clipped = stream_clipped.rio.clip([box(*dem_clipped.rio.bounds())], drop=True)

    # VALLEY DELINEATION
    flowline = flowline.simplify(10)
    points = get_cross_section_points(flowline, xs_spacing=200, xs_width=500, xs_point_spacing=10)
    points.index = range(len(points))

    points = add_attributes_to_xs(wbt, points, dem_clipped, flow_dir_clipped, hillslope_clipped)
    

   break_points = []
   for xs_id in points['cross_section_id'].unique():
       df = points.loc[points['cross_section_id'] == xs_id]
       df = rezero_xsection(df)

       # positive break point
       df_pos = df.loc[df['alpha'] >= 0]
       pos_break_point = find_break_point(df_pos)

       # negative break point
       df_neg = df.loc[df['alpha'] <= 0]
       neg_break_point = find_break_point(df_neg)
       break_points.append(pos_break_point)
       break_points.append(neg_break_point)
       # WRONG, returning index in xs not in points!
       
   break_points = [bp for bp in break_points if bp is not None]
   break_points = points.loc[break_points]

   valley_polygon = points_to_polygon(points, hillslope_clipped)

