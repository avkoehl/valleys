"""
input is clipped:
    conditioned_dem raster, 
    flow direction raster,
    stream raster, 
    stream flowline linestring
    hillslope raster

output:
    polygon
"""
import os

import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray
import scipy  
from scipy import signal
import xarray as xr

from valleys.cross_section import get_cross_section_points
from valleys.terrain import compute_terrain_rasters

def rezero_xsection(xs_points):
    # if the stream centerline was smoothed, 
    # then the starting point (alpha = 0) may not be where the stream was (hand = 0)
    # need to find that point, and adjust the alpha values accordingly 
    # alternatively could use the streamline and find the nearest point
    temp = xs_points.copy()

    # in theory, if xs points are very close they can be the same cell
    # this can cause the data to be the same, leading to a failure in peak finding
    # so we need to remove 'duplicate' points
    # better would be to drop duplicates based on the pointx and pointy but this is a good start
    # temp = temp.drop_duplicates(subset=['elevation', 'slope', 'curvature', 'hand'])
    # TODO: think on this, as doing this may remove alpha == 0

    if temp.loc[temp['alpha'] == 0]['hand'].iloc[0] != min(temp['hand']):
        min_ind = temp['hand'].idxmin(skipna=True)        
        temp['alpha'] = (temp['alpha'] - temp['alpha'][min_ind])

    return temp

def find_break_point(df, peak_threshold=0.002, slope_threshold=.4, hand_threshold=100):
    """ df is a half profile dataframe """

    # make sure that the alphas are in the right order (start at 0)
    df = df.reindex(df['alpha'].abs().sort_values().index)
    inds, _ =  signal.find_peaks(-df['curvature'], height=peak_threshold)

    if len(inds):
        key_point_inds = [0] + list(inds) + [len(df) - 1]
        key_points = df.iloc[key_point_inds]

        # hand threshold
        exceeds = (key_points['hand'] > hand_threshold).idxmax()
        alpha_exceeds = key_points['alpha'].loc[exceeds]
        key_points = key_points.loc[key_points['alpha'].abs() < np.abs(alpha_exceeds)]

        if not len(key_points):
            return (None, df['point_id'].iloc[inds].to_list())

        key_points['profile_step_slope'] = key_points['hand'].diff(periods=-1).abs() / key_points['alpha'].diff(periods=-1).abs()

        # find the first point where the slope is greater than the threshold
        exceeds = (key_points['profile_step_slope'] > slope_threshold)
        if exceeds.any():

            # if the first point is the channel point, then the break point is the point right after it
            if exceeds.iloc[0]:
                break_point = df['point_id'].iloc[exceeds.argmax() + 1]
            else:
                # return the first point that exceeds the threshold
                break_point = df['point_id'].iloc[exceeds.argmax()]
            return (break_point, df['point_id'].iloc[inds].to_list())

    # if nothing found return None
    return (None, df['point_id'].iloc[inds].to_list()) 

def determine_threshold(break_points):
    threshold = break_points['hand'].quantile(.6)
    return threshold

def rioxarray_sample_points(raster, points, method='nearest'):
    """ Sample points from raster using rioxarray """
    xs = xr.DataArray(points.geometry.x.values, dims='z')
    ys = xr.DataArray(points.geometry.y.values, dims='z')
    values = raster.sel(x=xs, y=ys, method=method).values
    return values

def add_attributes_to_xs(wbt, points, dem, stream, flow_dir, hillslope):
    # save rasters to temp files 
    dem_file = os.path.join(wbt.work_dir, 'temp_dem.tif')
    stream_file = os.path.join(wbt.work_dir, 'temp_stream.tif')
    flow_dir_file = os.path.join(wbt.work_dir, 'temp_flow_dir.tif')

    dem.rio.to_raster(dem_file)
    stream.rio.to_raster(stream_file)
    flow_dir.rio.to_raster(flow_dir_file)

    rasters = compute_terrain_rasters(wbt, dem_file, stream_file, flow_dir_file, sigma=1.2)
    dem = rioxarray.open_rasterio(dem_file, masked=True).squeeze()
    slope = rioxarray.open_rasterio(rasters['slope'], masked=True).squeeze()
    curvature = rioxarray.open_rasterio(rasters['curvature'], masked=True).squeeze()
    hand = rioxarray.open_rasterio(rasters['hand'], masked=True).squeeze()


    points['point_id'] = np.arange(len(points))

    points['elevation'] = rioxarray_sample_points(dem, points)
    points['slope'] = rioxarray_sample_points(slope, points)
    points['curvature'] = rioxarray_sample_points(curvature, points)
    points['hand'] = rioxarray_sample_points(hand, points)
    points['hillslope'] = rioxarray_sample_points(hillslope, points)

    points = points.loc[~points['elevation'].isna()]
    points = points.loc[~points['slope'].isna()]


    # cleanup temp files
    os.remove(dem_file)
    os.remove(stream_file)
    os.remove(flow_dir_file)
    return points


def polygonize(raster):
    raster.rio.to_raster('temp.tif', dtype=np.uint8)
    with rasterio.open('temp.tif') as src:
        raster_array = src.read(1)
        mask = raster_array == 1
        
        polygons = []
        for geom, value in rasterio.features.shapes(raster_array, mask=mask, transform=src.transform):
            if value == 1:  #
                polygons.append(shape(geom))
    os.remove('temp.tif')
    return polygons
    
def points_to_polygon(break_points, threshold, hand_raster, slope_raster, 
                      slope_threshold=.4, buffer=1):

    hand = hand.where(hand != -32768)
    values = hand.where(hand < (threshold + buffer))
    
    # binarize
    values = values.where(~np.isnan(values), 0)
    values = (values > 0).astype(int)
    
    # fill
    values.data = scipy.ndimage.binary_fill_holes(values.data)
    values = values.where(values != 0)
    
    # filter by slope
    values = values.where((slope < slope_threshold) & (values == 1))
    values.data = values.data.astype(np.uint8)

    # polygonize
    polygons = polygonize(values)
    polygons = gpd.GeoDataFrame(geometry=polygons, crs=3310)
    polygons['geometry'] = polygons['geometry'].apply(close_holes)
    return polygons

def close_holes(poly):
        if len(poly.interiors):
            return Polygon(list(poly.exterior.coords))
        return poly

def get_break_points(points):
    """
    Find the break points for each cross section
    Can return None, One, or Two break points

    cross_section_id, pos_break_point_id, neg_break_point_id, pos_peak_points, neg_peak_points

    """

    break_points = []
    for xs_id in points['cross_section_id'].unique():
        print(xs_id)
        df = points.loc[points['cross_section_id'] == xs_id]
        df_pos = df.loc[df['alpha'] >= 0]
        df_neg = df.loc[df['alpha'] <= 0]

        if len(df_pos) < 5 or len(df_neg) < 5:
            break_points.append((xs_id, None, None, None, None))
            continue

        df = rezero_xsection(df)

        pos_break_point, pos_inds = find_break_point(df_pos)
        neg_break_point, neg_inds = find_break_point(df_neg)

        break_points.append(
            (xs_id, pos_break_point, neg_break_point, pos_inds, neg_inds))
    break_points = pd.DataFrame(break_points, columns=['cross_section_id', 'pos_break_point_id', 'neg_break_point_id', 'pos_peak_points', 'neg_peak_points'])
    return break_points

def delineate_valley(wbt, dem_raster, flowline, stream_raster, hillslope_raster, flow_dir_raster):
    flowline = flowline.simplify(10)
    points = get_cross_section_points(flowline, xs_spacing=100, xs_width=200, xs_point_spacing=10)

    points = add_attributes_to_xs(wbt, points, dem_raster, stream_raster, flow_dir_raster)
    break_points = get_break_points(points)

    threshold = determine_threshold(break_points, hillslope)
    valley_polygon = points_to_polygon(points, threshold, hand, slope)

    return {'geometry': valley_polygon, 'hand_threshold': threshold}
