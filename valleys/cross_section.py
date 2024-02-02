"""
Code for methods that depend on cross section analysis

1. simplify stream linestring
2. get points on simplified stream linestring
3. get points on perpendicular lines to the simplified stream linestring
"""

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString, Point, Polygon
import rioxarray

def get_cross_section_points(linestring, simplify=True, tolerance=20, xs_spacing=5, 
                             xs_width=100, xs_point_spacing=10, crs=3310):
    """ Get Cross Section Points 
    """
    linestring = linestring.simplify(tolerance) if simplify else linestring

    # get points on the linestring
    points = _get_points_on_linestring(linestring, xs_spacing)
    
    alphas = list(range(-xs_width, xs_width+xs_point_spacing, xs_point_spacing))

    # for each point sample points on either side of the linestring
    dfs = []
    for i,point in enumerate(points):
        A, B = _get_nearest_vertices(point, linestring)
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
    """ Get Points on Linestring """
    points = []
    for i in range(0, int(linestring.length), spacing):
        point = linestring.interpolate(i)
        points.append(point)
    return points

    return points_df

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
