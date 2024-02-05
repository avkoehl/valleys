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
import rasterio
import rioxarray
import scipy  
from shapely.geometry import shape, MultiPolygon, Polygon
from scipy import signal
import xarray as xr
import matplotlib.pyplot as plt

from valleys.cross_section import get_cross_section_points
from valleys.terrain import compute_terrain_rasters

def filter_ridge_crossing(profile):
    # check if the profile crosses a ridgeline
    # if a point has positive curvature, and there is a big jump in HAND relative to the change in elevation
    # then it is likely that the profile has crossed a ridgeline
    # remove the points after the ridgeline
    # sometimes seems that curvature is not always positive when crossing a ridgeline
    temp = profile.copy()
    temp['diff'] = temp['alpha'].abs()
    # sort on diff column
    temp = temp.sort_values(by='diff', ascending=True)

    temp['hand_diff'] = temp['hand'].diff()
    temp['elev_diff'] = temp['elevation'].diff()
    temp['ratio'] = temp['hand_diff'] / temp['elev_diff']
    temp['ratio'] = temp['ratio'].abs()

    # find the first point where the ratio is greater than 2
    for i in range(len(temp)):
        if temp['ratio'].iloc[i] > 3 and temp['hand_diff'].iloc[i] > 20:
            #            if temp['curvature'].iloc[i] > 0:
            keep = temp.iloc[:i]
            return profile.loc[profile['point_id'].isin(keep['point_id'])]
    return profile

def prepare_xsection(xs_points):
    """
    This function prepares the cross section points for analysis
    0. make sure that the cross section has the correct columns
    1. remove any duplicate points
    2. make sure that there are enough points
    3. recenter the alpha values so that the channel is at 0
    4. make sure there are points on both sides of the channel
    5. remove points that are passed a ridgeline
    6. double check that there are enough points

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
    elif len(temp['strm_val'].dropna().unique()) == 1:
        strm_val = temp['strm_val'].dropna().unique()[0]
        min_point = temp['hand'].loc[temp['strm_val'] == strm_val].idxmin()
        min_alpha = temp['alpha'].loc[min_point]
        temp['alpha'] = temp['alpha'] - min_alpha
    elif len(temp['strm_val'].dropna().unique()) > 1:
        stream_points = temp.loc[temp['strm_val'].notna()]
        # find stream point with min abs alpha
        min_point = stream_points['alpha'].abs().idxmin()
        min_alpha = stream_points['alpha'].loc[min_point]
        temp['alpha'] = temp['alpha'] - min_alpha

    # 4. make sure there are points on both sides of the channel
    if len(temp.loc[temp['alpha'] < 0]) == 0 or len(temp.loc[temp['alpha'] > 0]) == 0:
        return None

    # 5. remove points that are passed a ridgeline
    temp_neg = filter_ridge_crossing(temp.loc[temp['alpha'] <= 0])
    temp_pos = filter_ridge_crossing(temp.loc[temp['alpha'] >= 0])
    temp = pd.concat([temp_neg, temp_pos])
    # remove duplicate points
    temp = temp.drop_duplicates(subset=['point_id'])

    # 6. double check that there are enough points
    if len(temp) < 5:
        return None

    if len(temp.loc[temp['alpha'] < 0]) == 0 or len(temp.loc[temp['alpha'] > 0]) == 0:
        return None
    return temp

def find_half_profile_break_point(profile, peak_ids, slope_threshold=20):
    profile = profile.copy()
    profile = profile.reindex(profile['alpha'].abs().sort_values().index)
    profile['int_index'] = np.arange(len(profile))

    # add the first and last points
    key_point_ids = peak_ids.copy()
    key_point_ids.insert(0, profile.iloc[0]['point_id'])
    key_point_ids.append(profile.iloc[-1]['point_id'])
    key_points = profile.loc[profile['point_id'].isin(key_point_ids)]
    key_points = key_points.drop_duplicates(subset=['alpha'])

    # get mean slope between each pair of points
    # foreach pair of points, get all the points between them and get mean of slope column
    slopes = []
    for ind, (label, row) in enumerate(key_points.iterrows()):
        if ind == len(key_points) - 1:
            break
        point_id_1 = row['point_id']
        point_id_2 = key_points.iloc[ind+1]['point_id']

        pos1 = profile.loc[profile['point_id'] == point_id_1]['int_index'].values[0]
        pos2 = profile.loc[profile['point_id'] == point_id_2]['int_index'].values[0]

        points = profile.iloc[pos1:pos2]
        mean_slope = points['slope'].mean()
        slopes.append((point_id_1, mean_slope))
    slopes = pd.DataFrame(slopes, columns=['point_id', 'mean_slope'])

    # get first index where slope > slope_threshold
    # else return None
    exceeds = slopes['mean_slope'] > slope_threshold
    if exceeds.any():
        return slopes.loc[exceeds].iloc[0]['point_id']
    return 

def find_xs_break_points(df, peak_threshold=0.002, slope_threshold=20):
    """ df is the points for a single cross section 
    has columns: alpha, curvature, hand, point_id, cross_section_id, hillslope, slope, strm_val

    # where these are point_ids
    # return (pos_break_point, neg_break_point, peak_ids)
    """

    # remove duplicates, recenter alpha, and make sure there are enough points
    df = prepare_xsection(df)

    # if the cross section is not valid, return None
    if df is None:
        return (None, None, None)

    # find the peaks in the curvature profile
    # if done after splitting, may not find the peak at the first point in the series (alpha == 0)
    # wheras if series is considered as a whole, the peak might be found, so do it before splitting
    peak_inds, _ =  signal.find_peaks(-df['curvature'], height=peak_threshold)
    peak_ids = df['point_id'].iloc[peak_inds].tolist()

    pos_bp = find_half_profile_break_point(df.loc[df['alpha'] >=0], peak_ids, slope_threshold)
    neg_bp = find_half_profile_break_point(df.loc[df['alpha'] <= 0], peak_ids, slope_threshold)

    return (pos_bp, neg_bp, peak_ids)

def determine_threshold(break_points):
    threshold = break_points['hand'].quantile(.6)
    return threshold

def rioxarray_sample_points(raster, points, method='nearest'):
    """ Sample points from raster using rioxarray """
    xs = xr.DataArray(points.geometry.x.values, dims='z')
    ys = xr.DataArray(points.geometry.y.values, dims='z')
    values = raster.sel(x=xs, y=ys, method=method).values
    return values

def derive_attributes(wbt, dem, stream, flow_dir, hillslope):
    # save rasters to temp files 
    dem_file = os.path.join(wbt.work_dir, 'temp_dem.tif')
    stream_file = os.path.join(wbt.work_dir, 'temp_stream.tif')
    flow_dir_file = os.path.join(wbt.work_dir, 'temp_flow_dir.tif')

    dem.rio.to_raster(dem_file)
    stream.rio.to_raster(stream_file)
    flow_dir.rio.to_raster(flow_dir_file)

    rasters = compute_terrain_rasters(wbt, dem_file, stream_file, flow_dir_file, sigma=1.2)

    # cleanup temp files
    os.remove(dem_file)
    os.remove(stream_file)
    os.remove(flow_dir_file)
    return rasters

def add_attributes_to_xs(points, rasters, hillslope, stream):

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
    
def delineate_floor(threshold, hand):

    values = hand.where(hand < threshold) 
    values = values.where(np.isnan(values), 1)
    values = values.where(~np.isnan(values), 0)
    values = (values > 0).astype(int)
    
    # fill
    values.data = scipy.ndimage.binary_fill_holes(values.data)
    values = values.where(values != 0)
    
    # TODO: redo this filter, for now ignored because slope raster isn't quite right -
    # its too small, has different dimensions than the hand raster
    # in future will calculate terrain rasters on the basin not the subbasin
    # store as a dataset rather than a bunch of files
    # filter by slope
    # values = values.where((slope < slope_threshold) & (values == 1))
    #values.data = values.data.astype(np.uint8)

    # polygonize
    polygons = polygonize(values)
    polygons = gpd.GeoDataFrame(geometry=polygons, crs=3310)
    polygons['geometry'] = polygons['geometry'].apply(close_holes)
    # convert to multipolygon or single polygon
    if len(polygons) > 1:
        polygon = MultiPolygon(polygons['geometry'].values)
        return polygon, values
    return polygons['geometry'].iloc[0], values

def close_holes(poly):
    if len(poly.interiors):
        return Polygon(list(poly.exterior.coords))
    return poly

def get_break_points(points):
    """
    Find the break points for each cross section
    Can return None, One, or Two break points per cross section

    Returns dataframe
        cross_section_id, pos_break_point_id, neg_break_point_id, [peak_point_ids]

    """

    break_points = []
    for xs_id in points['cross_section_id'].unique():
        pos, neg, peaks = find_xs_break_points(points.loc[points['cross_section_id'] == xs_id])
        break_points.append(
            (xs_id, pos, neg, peaks))

    break_points_df = pd.DataFrame(break_points, columns=['cross_section_id', 'pos', 'neg', 'peak_ids'])
    combined = break_points_df['pos'].dropna().to_list() + break_points_df['neg'].dropna().to_list()

    break_points = points.loc[points['point_id'].isin(combined)]
    # break_points.to_file('break_points.shp')
    return break_points


def plot_xs_break_points(df, pos, neg, peaks):
    #xs_id = 3
    #df = points.loc[points['cross_section_id'] == xs_id]
    #pos, neg, peaks = find_xs_break_points(df)
    #plot_xs_break_points(df, pos, neg, peaks)
    ##        # df is full profile of xs
    df = prepare_xsection(df)
    fig, ax = plt.subplots(4,1)

    ax[0].scatter(df['alpha'], df['elevation'])
    ax[0].scatter(df.loc[df['point_id'].isin(peaks), 'alpha'], df.loc[df['point_id'].isin(peaks), 'elevation'], c='b')
    ax[0].scatter(df.loc[df['point_id'] == pos, 'alpha'], df.loc[df['point_id'] == pos, 'elevation'], c='r')
    ax[0].scatter(df.loc[df['point_id'] == neg, 'alpha'], df.loc[df['point_id'] == neg, 'elevation'], c='g')

    ax[1].scatter(df['alpha'], df['hand'])
    ax[1].scatter(df.loc[df['point_id'].isin(peaks), 'alpha'], df.loc[df['point_id'].isin(peaks), 'hand'], c='b')
    ax[1].scatter(df.loc[df['point_id'] == pos, 'alpha'], df.loc[df['point_id'] == pos, 'hand'], c='r')
    ax[1].scatter(df.loc[df['point_id'] == neg, 'alpha'], df.loc[df['point_id'] == neg, 'hand'], c='g')

    ax[2].plot(df['alpha'], df['slope'])
    ax[2].scatter(df.loc[df['point_id'].isin(peaks), 'alpha'], df.loc[df['point_id'].isin(peaks), 'slope'], c='b')
    ax[2].scatter(df.loc[df['point_id'] == pos, 'alpha'], df.loc[df['point_id'] == pos, 'slope'], c='r')
    ax[2].scatter(df.loc[df['point_id'] == neg, 'alpha'], df.loc[df['point_id'] == neg, 'slope'], c='g')

    ax[3].plot(df['alpha'], df['curvature'])
    ax[3].scatter(df.loc[df['point_id'].isin(peaks), 'alpha'], df.loc[df['point_id'].isin(peaks), 'curvature'], c='b')
    ax[3].scatter(df.loc[df['point_id'] == pos, 'alpha'], df.loc[df['point_id'] == pos, 'curvature'], c='r')
    ax[3].scatter(df.loc[df['point_id'] == neg, 'alpha'], df.loc[df['point_id'] == neg, 'curvature'], c='g')
    plt.show()

def delineate_valley(wbt, dem_raster, flowline, stream_raster, hillslope_raster, flow_dir_raster):
    flowline = flowline.simplify(10)
    points = get_cross_section_points(flowline, xs_spacing=50, xs_width=200, xs_point_spacing=10)

    points = add_attributes_to_xs(wbt, points, dem_raster, stream_raster, flow_dir_raster)
    break_points = get_break_points(points)

    threshold = determine_threshold(break_points, hillslope)
    valley_polygon = points_to_polygon(points, threshold, hand, slope)

    return {'geometry': valley_polygon, 'hand_threshold': threshold}
