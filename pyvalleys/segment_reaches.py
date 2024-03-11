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
import pandas as pd
import numpy as np
from shapely.geometry import LineString, Point
import shapely
import rioxarray
import ruptures as rpt

from pyvalleys.centerline import get_centerline
from pyvalleys.cross_section import get_cross_section_lines
from pyvalleys.gis import rioxarray_sample_points

def get_pour_points(valley_polygon, flowline, flow_acc, spacing, pen):
    centerline = get_centerline(valley_polygon, flowline, flow_acc)

    # make sure flowline going downhill, make sure centerline going downhill
    # make sure centerline's last point matches the last point of the flowline
    # need to extend its last vertex to the last vertex of the flowline
    # make sure cross section code always draws cross section at end of the line

    # get breakpoints, including last point and first point
    # create the associated line strings
    # segment_id, start_point, end_point, geometry [linestring], mean_width, lenght

    # if flowline < 1km: 
    # somehow filter/combine these so that no segment is less than 1km

    lines = get_width_series(valley_polygon, centerline, spacing=spacing)

    bps = detect_segment_breakpoints(lines, pen=pen)

    ppts = bps.intersection(flowline)
    points_list = []
    inds = []
    for ind,point in ppts.items():
        if isinstance(point, shapely.geometry.multipoint.MultiPoint):
            # just keep first one
            point = point.geoms[0]
        points_list.append(point)
        inds.append(ind)
    ppts = gpd.GeoSeries(points_list, index=inds, crs = bps.crs)
    ppts = ppts.loc[ppts.type == "Point"] # some might not have intersection

    start_ind = -1 
    locs = []
    for bpid in ppts.index:
        slines = lines.loc[(lines['csid'] > start_ind) & (lines['csid'] <= bpid) ]
        start_ind = bpid
        locs.append(slines.index)

    mean_widths = []
    for l in locs:
        rows = lines.loc[l]
        mean_width = rows['widths'].mean()
        mean_widths.append(mean_width)


    segment_lengths = []
    for point in ppts:
        # get next point
        # cut linestring at those two points, get length


    segment_ids = list(range(len(segment_lengths)))
    ppts = gpd.GeoDataFrame({
        'segment_id': segment_ids,
        'geometry': ppts,
        'mean_width': mean_widths,
        'length': segment_lengths})

    # add point for the last point on the flowline
    # double check that point has higher flow accumulation than first point
    # add segment_length
    # width = NA, length = NA

    # filter
    # if point very close to end -> remove

    # if any segment is less than 500m



    return ppts

def get_width_series(valley_polygon, centerline, spacing):
    lines = get_cross_section_lines(centerline, xs_spacing=spacing, xs_width = 4000)
    clipped_lines = []
    for ind,row in lines.iterrows():
        ls = row['geometry']
        center_point = ls.centroid

        ls = gpd.GeoSeries(ls, crs=lines.crs)
        ls = ls.clip(valley_polygon)
        ls = ls.explode(index_parts=False)
        ls = ls.loc[ls.intersects(center_point.buffer(3))]

        if len(ls) == 1:
            clipped_lines.append(
                    {'csid': row['cross_section_id'], 'geometry': ls.geometry.iloc[0]})

    df = pd.DataFrame(clipped_lines)
    df = gpd.GeoDataFrame(df, geometry='geometry', crs=lines.crs)
    df['widths'] = df.length
    return df
    
def detect_segment_breakpoints(lines, pen=10):
    signal = lines['widths'].values
    algo = rpt.Pelt(model='rbf').fit(signal)
    result = algo.predict(pen=pen)
    # result's last element is the length of the series

    # if none return just the last point
    # otherwise return all the breakpoints and the last point
    if len(result) == 1:
        return lines.iloc[[-1]]
    else:
        result = result[0:-1]
        result.append(len(lines)-1)
        return lines.iloc[result]
