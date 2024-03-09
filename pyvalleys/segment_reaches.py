"""
input: 
    valley_floors
    flow_dir
    flow_acc
    flowpaths

output:
    reach_flowpaths 
    reach_subbasins
    reach_hillslopes

steps:
    1. get centerlines
    2. get continous width measurements
    3. add points where significant width change
    4. remove points if segment too small
    5. clip points to flowpaths somehow
    6. watershed algorithm
"""

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, Point
import rioxarray
from shapelysmooth import taubin_smooth # prefer taubin unless need to preserve nodes
from shapelysmooth import chaikin_smooth

from pyvalleys.centerline import get_centerline
from pyvalleys.cross_section import get_cross_section_points
import ruptures as rpt

def get_width_series(valley_polygon, centerline, spacing=10):
    cs_points = get_cross_section_points(centerline, simplify=False, xs_spacing=spacing, xs_width = 2000, xs_point_spacing=2000)

    lines = []
    for csid in cs_points['cross_section_id'].unique():
        points = cs_points.loc[cs_points['cross_section_id'] == csid]
        start = points.iloc[0]
        end = points.iloc[-1]
        linestring = gpd.GeoSeries(LineString([start.point, end.point]), crs=3310)
        linestring = linestring.clip(valley_polygon)
        linestring = linestring.explode(index_parts=False)
        linestring = linestring.loc[linestring.intersects(centerline)]
        if len(linestring) == 1:
            lines.append(linestring.iloc[0])
    
    lines = gpd.GeoDataFrame(geometry=lines, crs=3310)
    lines['widths'] = lines.length
    return lines

def detect_segment_breakpoints(lines, pen=10):
    signal = lines['widths'].values
    algo = rpt.Pelt(model='rbf').fit(signal)
    result = algo.predict(pen=pen)
    return lines.iloc[result[0:-1]]

def locate_stream_coordinates(breakpoints, flowpath):
    # find where breakpoint line intersects flowpath linestring
    # if no intersect give error
    # if multiple intersects keep closest

    return coordinates

# try snapping coordinates to flow acc
# run watershed tool? in qgis with some sample points


base = "../testing/"
valleys = gpd.read_file(f"{base}valley_floors.shp")
flowlines = gpd.read_file(f"{base}/1805000203/terrain_attributes/flowpaths.shp")
flow_acc = rioxarray.open_rasterio(f"{base}/1805000203/terrain_attributes/flow_acc.tif", masked=True).squeeze()

sid = 9 #9
valley_polygon = valleys.loc[valleys['ID'] == sid]['geometry'].buffer(1).iloc[0]
flowline = flowlines.loc[flowlines['STRM_VAL'] == sid]['geometry'].iloc[0]

centerline = get_centerline(valley_polygon, flowline, flow_acc)
widths = get_width_series(valley_polygon, centerline)
bps = detect_segment_breakpoints(widths)

    




fig, ax = plt.subplots()
widths.plot(ax=ax)
bps.plot(ax=ax, color='red')
gpd.GeoSeries(valley_polygon).plot(ax=ax, facecolor='none', edgecolor='black')
gpd.GeoSeries(flowline).plot(ax=ax)
gpd.GeoSeries(centerline).plot(ax=ax, color='green', linestyle='-')
plt.show()
