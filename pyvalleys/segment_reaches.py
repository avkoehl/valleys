"""
input: 
    valley_floors
    flow_dir
    flow_acc
    flowpaths

output:
    reach_pour_points

steps:
    1. get centerlines
    2. get continous width measurements
    3. detect segments with mean shift of width series
    4. get points of the breakpoints clipped to flowline
"""

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, Point
import rioxarray

from pyvalleys.centerline import get_centerline
from pyvalleys.cross_section import get_cross_section_points
import ruptures as rpt

def get_pour_points(valley_polygon, flowline, flow_acc, spacing=10, pen=10):
    centerline = get_centerline(valley_polygon, flowline, flow_acc)

    lines = get_width_series(valley_polygon, centerline, spacing=spacing)
    bps = detect_segment_breakpoints(lines, pen=pen)
    ppts = bps.intersection(flowline)
    ppts = ppts.loc[ppts.type == "Point"]
    return ppts

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
    return lines.iloc[result] # this includes the last point in the series
