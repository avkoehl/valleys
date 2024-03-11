"""
Code for methods that depend on cross section analysis

1. get points on stream linestring at xs_spacing
2. get points on perpendicular lines to the stream linestring 
   at xs_point_spacing up to xs_width on either side of the stream
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely
from shapely.geometry import LineString, Point, MultiPoint
from shapely.ops import nearest_points

def get_cross_section_lines(linestring, xs_spacing, xs_width, crs=3310):
    points = _get_points_on_linestring(linestring, xs_spacing)
    width = int(xs_width / 2)

    dfs = []
    for i,point in enumerate(points):
        A, B = _get_nearest_vertices2(point, linestring)
        p1 = _sample_point_on_perpendicular_line(point, A, B, width)
        p2 = _sample_point_on_perpendicular_line(point, A, B, -width)
        xs = LineString([p1,p2])
        df = pd.DataFrame(
                [{'cross_section_id': i,
                 'geometry': xs}])
        dfs.append(df)
    df = pd.concat(dfs)
    lines = gpd.GeoDataFrame(df, geometry='geometry', crs=crs)
    lines = lines.reset_index(drop=True)
    return lines

def get_cross_section_points(linestring, xs_spacing=5, 
                             xs_width=100, xs_point_spacing=10, crs=3310):
    """ Get Cross Section Points 
    """
    # get points on the linestring
    points = _get_points_on_linestring(linestring, xs_spacing)
    
    alphas = list(range(-xs_width, xs_width+xs_point_spacing, xs_point_spacing))

    # for each point sample points on either side of the linestring
    dfs = []
    for i,point in enumerate(points):
        A, B = _get_nearest_vertices2(point, linestring)
        df = pd.DataFrame(
                {'alpha': alphas,
                 'point': [_sample_point_on_perpendicular_line(point, A, B, alpha) for alpha in alphas],
                 'cross_section_id': i})
        dfs.append(df)

    points_df = pd.concat(dfs)
    points_df = gpd.GeoDataFrame(points_df, geometry='point', crs=crs)
    return points_df

# --- Internal Functions ----------------------------------------------------- #
def _get_points_on_linestring(linestring, spacing):
    points = []
    for i in range(0, int(linestring.length), spacing):
        point = linestring.interpolate(i)
        points.append(point)
    return points

def _get_nearest_vertices2(point, linestring):
    # its known that the point is on the linestring
    linestring_copy = linestring
    linestring = gpd.GeoSeries(linestring).clip(point.buffer(5)).iloc[0]

    if (isinstance(linestring, shapely.geometry.multilinestring.MultiLineString)):
        ls = gpd.GeoSeries(linestring)
        ls = ls.explode(index_parts=False).reset_index(drop=True)
        closest = ls.loc[ls.distance(point).idxmin()]
        linestring = closest


    mp = MultiPoint(linestring.coords)
    nearest_point = nearest_points(point, mp)[1]
    index = list(linestring.coords).index((nearest_point.x, nearest_point.y))

    if index == 0:
        return [nearest_point, Point(linestring.coords[1])]

    if index == (len(linestring.coords) -1):
        return [nearest_point, Point(linestring.coords[-2])]

    second_nearest = min(Point(linestring.coords[index - 1]),
                         Point(linestring.coords[index + 1]),
                         key=point.distance)
    return [nearest_point, second_nearest]

def _get_nearest_vertices(point, linestring):
    line_coords = linestring.coords
    distances = [point.distance(Point(x, y)) for x, y in line_coords]
    nearest_indices = np.argsort(distances)[:2]
    nearest_vertices = [Point(line_coords[int(i)]) for i in nearest_indices]
    return nearest_vertices

def _sample_point_on_perpendicular_line(point, A, B, alpha):
    length = A.distance(B)
    x = point.x + alpha * (A.y - B.y) / length
    y = point.y + alpha * (B.x - A.x) / length
    return Point(x, y)
