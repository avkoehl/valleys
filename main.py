import os
import shutil

import geopandas as gpd
import rioxarray
from shapely.geometry import box
import pandas as pd

from valleys.utils import setup_wbt
from valleys.utils import chomp_raster
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
subbasins = rioxarray.open_rasterio(subbasin_files['subbasins'], masked=True).squeeze()
flowlines = gpd.read_file(stream_files['streams_shp'])
dem = rioxarray.open_rasterio(flow_acc_files['conditioned_dem'], masked=True).squeeze()
hillslopes = rioxarray.open_rasterio(subbasin_files['hillslopes'], masked=True).squeeze()
flow_dir = rioxarray.open_rasterio(flow_acc_files['flow_dir'], masked=True).squeeze()
stream = rioxarray.open_rasterio(stream_files['streams'], masked=True).squeeze()

# ------------ VALLEY DELINEATION ------------
# subbasin 1 xs 39 if spacing is 50 is example of cross section 

# work with sample
subbasin_id = 1

flowline = flowlines.loc[flowlines['STRM_VAL'] == subbasin_id]['geometry'].iloc[0]
dem, hillslope, flow_dir, stream = clip_to_subbasin(dem, hillslopes, flow_dir, stream, subbasins, subbasin_id)

points = get_cross_section_points(flowline, xs_spacing=50, xs_width=500, xs_point_spacing=10, tolerance=20)
points = add_attributes_to_xs(wbt, points, dem, stream, flow_dir, hillslope)

break_points_df = get_break_points(points)

# plot example
xs_id = 39
break_points = break_points_df.loc[break_points_df['cross_section_id'] == xs_id]


#---

bp = points.loc[points['point_id'] == break_points['pos_break_point_id'].iloc[0]]
peaks = break_points['pos_peak_points'].iloc[0] +  break_points['neg_peak_points'].iloc[0]
peaks = points.loc[points['point_id'].isin(peaks)]

fig, axes = plt.subplots(1,2)

x = points.loc[points['cross_section_id'] == xs_id]['alpha']
y = points.loc[points['cross_section_id'] == xs_id]['hand']

# add peaks
axes[0].plot(x,y)
axes[0].scatter(peaks['alpha'], peaks['hand'], c='r', label='peaks')
axes[0].scatter(bp['alpha'], bp['hand'], c='green', label='peaks')
axes[1].plot(x, points.loc[points['cross_section_id'] == xs_id]['curvature'])
axes[1].scatter(peaks['alpha'], peaks['curvature'], c='r', label='peaks')
plt.show()
