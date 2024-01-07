"""
Get features for the cells that are then sampled at the cross sections
"""
import os

import rasterio
import rioxarray
import geopandas as gpd


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
            min_ind = df['elevation'].idxmin()        
            temp.loc[temp['cross_section_id'] == ind, 'alpha'] = (df['alpha'] - df['alpha'][min_ind])

    return temp

def rioxarray_sample_points(raster, points, method='nearest'):
    """ Sample points from raster using rioxarray """
    xs = xr.DataArray(points.geometry.x.values, dims='z')
    ys = xr.DataArray(points.geometry.y.values, dims='z')
    values = raster.sel(x=xs, y=ys, method=method).values
    return values

def compute_terrain_features(wbt, dem_raster, sigma=1.5, log=False):
    # gaussian filter
    # slope
    # profile curvature
    # hand
    return

def get_terrain_features(points, dem, slope, curvature, hand, hillshade):
    # get the terrain features for the points
    # points should already be in the same crs as the rasters
    # points should have a 'cross_section_id' column
    # dem, slope, curvature, hand should be rioxarray objects

    points['elevation'] = rioxarray_sample_points(dem, points)
    points['slope'] = rioxarray_sample_points(slope, points)
    points['curvature'] = rioxarray_sample_points(curvature, points)
    points['hand'] = rioxarray_sample_points(hand, points)
    points['hillshade_id'] = rioxarray_sample_points(hillshade, points)
    points = rezero_alphas(points)

    return points
