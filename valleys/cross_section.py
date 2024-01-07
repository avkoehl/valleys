"""
Code for methods that depend on cross section analysis

1. get stream linestring
2. simplify stream linestring
3. get points on simplified stream linestring
4. Foreach point get points on either side of the stream linestring
    df: cross_section_id, alpha, point
5. for those points get elevation, slope, curvature
6. rezero alphas so that the lowest point (should be close to 0 - where the original stream is intersected is at 0)
5. save dataframe  cross_section_id, corrected_alpha, point, elevation, slope, curvature
"""

import os

from label_centerlines import get_centerline
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import LineString, Point, Polygon
import rioxarray


def vectorize_stream(wbt, stream_binary_raster):
    """ stream raster should contain a single stream """
    thinned = morphology.thin(stream_binary_raster.data)
    thinned = thinned.astype(int)

    with rasterio.open('temp.tif', 'w', driver='GTiff', height = thinned.shape[0],
                       width = thinned.shape[1], count=1, dtype=str(thinned.dtype),
                       crs=hillslope_raster.rio.crs, transform=hillslope_raster.rio.transform()) as dst:
        dst.write(thinned, 1)

    wbt.raster_to_vector_lines(os.path.abspath('temp.tif'), 'temp.shp')
    stream = gpd.read_file(os.path.join(wbt.work_dir, 'temp.shp'), crs=stream_binary_raster.rio.crs)

    os.remove('temp.tif')
    return

def generate_cross_section_lines(points):
    #TODO
    pass

def preprocess_channel(linestring, method="centerline", threshold=None, 
                       hand=None, contour_levels=None):
    """ Preprocess Channel 

    Its important that the channel meanders are minimized before
    extracting cross sections so that the cross sections are more
    likely to be orthogonal to the valley walls.

    Using the minimum_rotated_rectangle method is not ideal because 
    it will not work well with bends in the river around the hills - i.e
    the valley itself is sinous.

    """

    if method not in ["centerline", "simplify"]:
        raise ValueError("Invalid method")

    if method == "centerline":
        return _centerline(hand, linestring, contour_levels)

    if method == "simplify":
        return linestring.simplify(threshold)

def _centerline(hand, linestring, contour_levels=[0, 5, 10, 100, 150]):
    """
    Get centerline of the polygon representing the channel floodplain

    This polygon is extracted from the contour level of the DEM
    using matplotlib
    """
    extent = _get_extent(hand)
    cs = plt.contourf(hand, levels=contour_levels, colors='black', 
                    extent=extent, origin='upper')

    # get polygon
    polygon = None
    for collection in cs.collections:
        if polygon is not None:
            break
        for path in collection.get_paths():
            poly = Polygon(path.vertices)
            if poly.intersects(linestring):
                polygon = poly
                break

    centerline = get_centerline(polygon)
    # consider extending the centerline to the edges of the polygon

    return centerline

def _get_points_on_linestring(linestring, spacing):
    """ Get Points on Linestring """
    points = []
    for i in range(0, int(linestring.length), spacing):
        point = linestring.interpolate(i)
        points.append(point)
    return points

def get_cross_section_points(linestring, xs_spacing=5, 
                             xs_width=100, xs_point_spacing=10, crs=3310):
    """ Get Cross Section Points """

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
