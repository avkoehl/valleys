"""
Get features for the cells that are then sampled at the cross sections
"""
import os

import geopandas as gpd
import numpy as np
import rasterio
import rioxarray
import xarray as xr

from valleys.utils import get_stream_raster_from_hillslope


def rezero_alphas(points):
    # if the stream centerline was smoothed, 
    # then the starting point (alpha == 0) may not be where the stream was (elevation = 0)
    # need to find that point, and adjust the alpha values accordingly 

    # could use the streamline and find the nearest point
    temp = points.copy()

    offsets = {}
    for ind in temp['cross_section_id'].unique():
        df = temp.loc[points['cross_section_id'] == ind]

        if df.loc[df['alpha'] == 0]['elevation'].iloc[0] != min(df['elevation']):
            min_ind = df['elevation'].idxmin(skipna=True)        
            temp.loc[temp['cross_section_id'] == ind, 'alpha'] = (df['alpha'] - df['alpha'][min_ind])

    return temp


def rioxarray_sample_points(raster, points, method='nearest'):
    """ Sample points from raster using rioxarray """
    xs = xr.DataArray(points.geometry.x.values, dims='z')
    ys = xr.DataArray(points.geometry.y.values, dims='z')
    values = raster.sel(x=xs, y=ys, method=method).values
    return values

def compute_hand(wbt, dem_raster, hillslope):
    stream = get_stream_raster_from_hillslope(wbt, hillslope)
    stream.rio.to_raster(os.path.join(wbt.work_dir, 'stream.tif'))
    dem_raster.rio.to_raster(os.path.join(wbt.work_dir, 'dem.tif'))

    # get stream from hillslope for hand
    # need stream.tif to be 1 for stream, 0 for not stream
    # need to fill dem first
    wbt.fill_depressions('dem.tif', 'filled_dem.tif')
    wbt.elevation_above_stream('filled_dem.tif', 'stream.tif', 'hand.tif')
    hand = rioxarray.open_rasterio(os.path.join(wbt.work_dir, 'hand.tif'))
    hand = hand.squeeze()
    return hand

def _compute_terrain_features(wbt, dem_raster, hillslope, sigma=1.5, log=False):
    dem_raster.rio.to_raster(os.path.join(wbt.work_dir, 'dem.tif'))
    hand = compute_hand(wbt, dem_raster, hillslope)

    # gaussian filter
    wbt.gaussian_filter('dem.tif', 'smoothed.tif', sigma)
    wbt.slope('smoothed.tif', 'slope.tif')
    wbt.profile_curvature('smoothed.tif', 'curvature.tif')

    slope = rioxarray.open_rasterio(os.path.join(wbt.work_dir, 'slope.tif'))
    slope = slope.squeeze()
    slope = slope.where(slope != -32768)

    curvature = rioxarray.open_rasterio(os.path.join(wbt.work_dir, 'curvature.tif'))
    curvature = curvature.squeeze()
    curvature = curvature.where(curvature != -32768)
    return slope, curvature, hand

def get_terrain_features(wbt, points, dem, hillslope):
    # get the terrain features for the points
    # points should already be in the same crs as the rasters
    # points should have a 'cross_section_id' column
    # dem, slope, curvature, hand should be rioxarray objects

    slope, curvature, hand = _compute_terrain_features(wbt, dem, hillslope, sigma=1.5, log=False)

    points['elevation'] = rioxarray_sample_points(dem, points)
    points['slope'] = rioxarray_sample_points(slope, points)
    points['curvature'] = rioxarray_sample_points(curvature, points)
    points['hand'] = rioxarray_sample_points(hand, points)
    points['hillslope_id'] = rioxarray_sample_points(hillslope, points)
    points = rezero_alphas(points)

    # remove points where any value is nan or -32768
    points = points.replace(-32768, np.nan)
    points = points.dropna()

    # add unique id for each point
    points['point_id'] = np.arange(len(points))

    return points, slope, curvature, hand
