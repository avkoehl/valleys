"""
Code for methods that depend on cross section analysis
Methods that use regions/buffers/contours are in: region.py

Methods:
    - Plotting
        - plot_cross_section_profile     ✓
        - map_elevation_contours         ✓
        - map_cross_sections             ✓

    - Cross Section
        - preprocess_channel             ✓
        - get_cross_section_points       ✓


"""
from label_centerlines import get_centerline
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import shapely
from shapely.geometry import LineString, Point, Polygon
import rioxarray
import xrspatial

# ---------------- PLOTTING ---------------- #
def plot_cross_section_profile(points_df, ax=None):
    """ Plot Cross Section Profile """
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    df = points_df.copy()
    df = df.loc[np.isfinite(df['elevation'])]
    ax.scatter(df['alpha'], df['elevation'])
    ax.scatter(df['alpha'], df['slope'])
    return ax

def _generate_cross_section_lines(points_df):
    """ Generate Cross Section Lines """
    lines = []
    for index in points_df['cross_section'].unique():
        df = points_df.loc[points_df['cross_section'] == index]
        df = df.loc[np.isfinite(df['elevation'])]
        min_alpha_index = np.argmin(df['alpha'])
        max_alpha_index = np.argmax(df['alpha'])
        start = df['point'].iloc[min_alpha_index]
        end = df['point'].iloc[max_alpha_index]
        line = LineString([start, end])
        lines.append(line)
    lines = gpd.GeoDataFrame(geometry=lines)
    return lines

def _get_extent(raster):
    # need tuple of (left, right, bottom, top)
    # raster bounds is in format (left, bottom, right, top)
    bounds = raster.rio.bounds()
    return (bounds[0], bounds[2], bounds[1], bounds[3])

def map_elevation_contours(dem, ax=None, levels=20):
    """ Plot DEM, usually HAND with Contours on Top """
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    cs = ax.contourf(hand, cmap='terrain', levels=levels,
                     extent=_get_extent(hand), origin='upper')
    ax.clabel(cs, inline=True, fontsize=10)
    return ax


def map_cross_sections(points_df, hand, linestring, ax=None):
    """ Plot DEM with Cross Sections on Top """
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    cs = ax.contourf(hand, cmap='terrain', levels=20,
                     extent=_get_extent(hand), origin='upper')
    ax.clabel(cs, inline=True, fontsize=10)

    lines = _generate_cross_section_lines(points_df)
    lines.plot(ax=ax, color='orange')

    if not isinstance(linestring, gpd.GeoDataFrame):
        linestring = gpd.GeoDataFrame(geometry=[linestring])
    linestring.plot(ax=ax, color='red')

    return ax

# ---------------- CROSS SECTION ---------------- #

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
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    extent = _get_extent(hand)
    cs = ax.contourf(hand, levels=contour_levels, colors='black', 
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

def get_cross_section_points(linestring, hand, slope, xs_spacing=5, 
                             xs_width=100, xs_point_spacing=10):
    """ Get Cross Section Points """

    # get points on the linestring
    points = []
    for i in range(0, int(linestring.length), xs_spacing):
        point = linestring.interpolate(i)
        points.append(point)

    # for each point sample points on either side of the linestring
    dfs = []
    for i,point in enumerate(points):
        res = {}
        alphas, elevations, points = _sample_points(point, linestring, xs_width, xs_point_spacing, hand)
        _, slopes, _ = _sample_points(point, linestring, xs_width, xs_point_spacing, slope)
        res['alpha'] = alphas
        res['elevation'] = elevations
        res['slope'] = slopes
        res['point'] = points
        res['cross_section'] = i
        df = pd.DataFrame(res)
        dfs.append(df)

    points_df = pd.concat(dfs)

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

def _sample_points(point, linestring, width, spacing, raster):
    A,B = _get_nearest_vertices(point, linestring)
    npoints = width // spacing
    positive_alphas = [i * spacing for i in range(npoints+1)]
    negative_alphas = [-i for i in positive_alphas[::-1]]
    alphas = negative_alphas + positive_alphas
    points = [_sample_point_on_perpendicular_line(point, A, B, alpha) for alpha in alphas]
    values = [raster.sel(x=point.x, y=point.y, method='nearest').values.item() for point in points]
    return (alphas, values, points)
