"""
# slope
# curvature
# curvature of profile
# wavelet transform

# breakpoints to polygon
    # breakpoints filter
    # contour of hand with buffer
    # remove cells < slope theshold
    # fill holes
    # vectorize
"""
import os

import geopandas as gpd
import numpy as np
import rasterio
import scipy
from scipy import signal
from shapely.geometry import shape, Polygon


def get_half_profiles(df, xs_id):
    """ df is the points dataframe """
    df = df.loc[df['cross_section_id'] == xs_id]

    # negative alphas
    neg = df.loc[df['alpha'] <= 0]
    neg.loc[:, 'alpha'] = neg['alpha'].abs()
    neg = neg.sort_values(by='alpha', ascending=True)

    # positive alphas
    pos = df.loc[df['alpha'] >= 0]
    pos = pos.sort_values(by='alpha', ascending=True)

    return neg, pos

def best_peak_points(inds, df, slope_threshold):
    """ df is a half profile dataframe """
    for i,ind in enumerate(inds):
        if i == len(inds) - 1:
           slope = (df['elevation'].iloc[-1] - df['elevation'].iloc[ind]) / (df['alpha'].iloc[-1] - df['alpha'].iloc[ind])
        else:
            slope = (df['elevation'].iloc[inds[i+1]] - df['elevation'].iloc[ind]) / (df['alpha'].iloc[inds[i+1]] - df['alpha'].iloc[ind])

        if slope >= slope_threshold:
            return ind
    return 

def get_break_point(df, peak_threshold=0.005, slope_threshold=.25):
    """ df is a half profile dataframe """
    peak_inds, _ =  signal.find_peaks(-df['curvature'], height=peak_threshold)
    best_peak_ind = best_peak_points(peak_inds, df, slope_threshold)
    if best_peak_ind is not None:
        return df['point_id'].iloc[best_peak_ind]
    return

def get_break_points(df, peak_threshold=0.005, slope_threshold=.25):
    """ df is dataframe with all the points """
    point_ids = []
    for xs_id in df['cross_section_id'].unique():
        neg, pos = get_half_profiles(df, xs_id)
        neg_break = get_break_point(neg, peak_threshold, slope_threshold)
        pos_break = get_break_point(pos, peak_threshold, slope_threshold)
        point_ids.append(neg_break)
        point_ids.append(pos_break)

    # remove None values
    point_ids = [x for x in point_ids if x is not None]
    break_points = df.loc[df['point_id'].isin(point_ids)]
    return break_points

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
    
def map_valley_floor(break_points, hand, slope, quantile=.80, buffer=2, slope_threshold=25):

    # get contour
    #  TODO: remove outliers
    heights = break_points['hand']
    contour_height = heights.quantile(quantile)
    hand = hand.where(hand != -32768)
    values = hand.where(hand < (contour_height + buffer))
    
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

