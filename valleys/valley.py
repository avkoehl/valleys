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

def prepare_xsection(xs_points):
    """
    This function prepares the cross section points for analysis
    0. make sure that the cross section has the correct columns
    1. remove any duplicate points
    2. make sure that there are enough points
    3. recenter the alpha values so that the channel is at 0

    The cross section points are sampled along a line perpendicular to the 
    centerline, which is generally a smoothed version of the stream centerline
    as a result, the centerpoint of the cross section may not fall on the stream

    additionally, the cross section may not have enough points to be accurately
    useful in the peak detection analysis

    the cross section may also have duplicated points if for example the point spacing is smaller than
    the dem resolution and method='nearest' was used for the sampling

    recenter alpha so that the channel is at 0
      this is tricky because the cross section may intersect with the channel at several points
      alternatively, if the spacing of points is large, then the channel may be skipped entirely
      in which case we want to recenter to the lowest point
      in either case we want to recenter to the lowest point that ALSO is closest to the middle of the 
      cross section. Consider a case where the cross section is so wide that it goes into another valley
    """
    temp = xs_points.copy()

    # 0. make sure that the cross section has the correct columns
    req = ['alpha', 'curvature', 'hand', 'point_id', 'cross_section_id', 'hillslope', 'slope', 'strm_val']
    if not all([col in xs_points.columns for col in req]):
        raise ValueError(f'Columns missing from cross section points. Required: {req}')

    # 1. remove any duplicate points
    #  duplicated points can occur if the point spacing is smaller than the dem resolution
    #  alpha, point_id, and point will all be unique
    duplicted = temp.duplicated(subset=['elevation', 'slope', 'curvature', 'hand', 'hillslope', 'strm_val'])
    temp = temp.loc[~duplicted]

    # 2. make sure that there are enough points
    if len(temp) < 5:
        return None

    # 3. recenter the alpha values so that the channel is at 0
    #  Aprroach
    #   Check if multiple stream values are present in the cross section
    #          If yes, something wrong with the stream clipping step

    #   Case 1: No point has the stream value
    #        Recenter the alpha values to the point with the min(hand) value
    #        This isn't ideal, but its a good start
    #   Case 2: Only one point has the stream value
    #        Recenter the alpha values to that point
    #   Case 3: Multiple points have the stream value
    #         Find the point with the min(abs(alpha)) and recenter the alpha values to that point
    
    if len(temp['strm_val'].dropna().unique()) > 1:
        raise ValueError('Cross section has multiple stream values')

    if temp['strm_val'].isna().all():
        min_point = temp['hand'].idxmin()
        min_alpha = temp['alpha'].loc[min_point]
        temp['alpha'] = temp['alpha'] - min_alpha
    else if len(temp['strm_val'].dropna().unique()) == 1:
        strm_val = temp['strm_val'].dropna().unique()[0]
        min_point = temp['hand'].loc[temp['strm_val'] == strm_val].idxmin()
        min_alpha = temp['alpha'].loc[min_point]
        temp['alpha'] = temp['alpha'] - min_alpha
    else if len(temp['strm_val'].dropna().unique()) > 1:
        stream_points = temp.loc[temp['strm_val'].notna()]
        # find stream point with min abs alpha
        min_point = stream_points['alpha'].abs().idxmin()
        min_alpha = stream_points['alpha'].loc[min_point]
        temp['alpha'] = temp['alpha'] - min_alpha

    return temp

def find_half_profile_break_point(profile, peak_ids, slope_threshold=20):
    profile = profile.reindex(profile['alpha'].abs().sort_values().index)
    key_points = profile.loc[profile['point_id'].isin(peak_ids)]

    # add the first and last points
    key_points = key_points.append(profile.iloc[0])
    key_points = key_points.append(profile.iloc[-1])
    # remove duplicates
    key_points = key_points.drop_duplicates(subset=['alpha'])

    # get mean slope between each pair of points
    # foreach pair of points, get all the points between them and get mean of slope column
    slopes = []
    for ind, (label, row) in enumerate(key_points.iterrows()):
        point1 = row
        point2 = key_points.iloc[ind+1]

        # positions in the profile
        pos1 = profile.index.get_loc(point1['point_id'])
        pos2 = profile.index.get_loc(point2['point_id'])
        points = profile.iloc[pos1:pos2]
        mean_slope = points['slope'].mean()
        slopes.append(point1['point_id'], mean_slope)
    slopes = pd.DataFrame(slopes, columns=['point_id', 'mean_slope'])

    # get first index where slope > slope_threshold
    # else return None
    exceeds = slopes['mean_slope'] > slope_threshold
    if exceeds.any():
        return slopes.loc[exceeds].iloc[0]
    return 

def find_xs_break_points(df, peak_threshold=0.002, slope_threshold=20)
    """ df is the points for a single cross section 
    has columns: alpha, curvature, hand, point_id, cross_section_id, hillslope, slope, strm_val

    # where these are point_ids
    # return (pos_break_point, neg_break_point, peak_ids)
    """

    # remove duplicates, recenter alpha, and make sure there are enough points
    df = prepare_xsection(df)
    df['dist'] = df['alpha'].abs()

    # if the cross section is not valid, return None
    if df is None:
        return (None, None)

    # find the peaks in the curvature profile
    # if done after splitting, may not find the peak at the first point in the series (alpha == 0)
    # wheras if series is considered as a whole, the peak might be found, so do it before splitting
    peak_inds, _ =  signal.find_peaks(-df['curvature'], height=peak_threshold)
    peak_ids = df['point_id'].iloc[peak_inds]

    pos_bp = find_half_profile_break_point(df.loc[df['alpha'] >=0], peak_ids, slope_threshold)
    neg_bp = find_half_profile_break_point(df.loc[df['alpha'] <= 0], peak_ids, slope_threshold)

    return (pos_bp, neg_bp, peak_ids.tolist())

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
    dem = rioxarray.open_rasterio(rasters['dem_gauss'], masked=True).squeeze()
    slope = rioxarray.open_rasterio(rasters['slope'], masked=True).squeeze()
    curvature = rioxarray.open_rasterio(rasters['curvature'], masked=True).squeeze()
    hand = rioxarray.open_rasterio(rasters['hand'], masked=True).squeeze()


    points['point_id'] = np.arange(len(points))

    points['elevation'] = rioxarray_sample_points(dem, points)
    points['slope'] = rioxarray_sample_points(slope, points)
    points['curvature'] = rioxarray_sample_points(curvature, points)
    points['hand'] = rioxarray_sample_points(hand, points)
    points['hillslope'] = rioxarray_sample_points(hillslope, points)
    points['strm_val'] = rioxarray_sample_points(stream, points)

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

# TODO: redo to use the modified find_xs_break_points function
def get_break_points(points):
    """
    Find the break points for each cross section
    Can return None, One, or Two break points per cross section

    Returns dataframe
        cross_section_id, pos_break_point_id, neg_break_point_id, [peak_point_ids]

    """

    break_points = []
    for xs_id in points['cross_section_id'].unique():
        pos, neg, peaks = find_xs_break_points(points, xs_id)
        break_points.append(
            (xs_id, pos, neg, peaks)

    break_points = pd.DataFrame(break_points, columns=['cross_section_id', 'pos', 'neg', 'peak_ids'])
    return break_points

def delineate_valley(wbt, dem_raster, flowline, stream_raster, hillslope_raster, flow_dir_raster):
    flowline = flowline.simplify(10)
    points = get_cross_section_points(flowline, xs_spacing=50, xs_width=200, xs_point_spacing=10)

    points = add_attributes_to_xs(wbt, points, dem_raster, stream_raster, flow_dir_raster)
    break_points = get_break_points(points)

    threshold = determine_threshold(break_points, hillslope)
    valley_polygon = points_to_polygon(points, threshold, hand, slope)

    return {'geometry': valley_polygon, 'hand_threshold': threshold}
