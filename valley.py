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
import pandas as pd
import rioxarray
import scipy  
from scipy import signal
import xarray as xr

from valleys.cross_section import get_cross_section_points
from valleys.fun import compute_terrain_rasters

def rezero_xsection(xs_points):
    # if the stream centerline was smoothed, 
    # then the starting point (alpha == 0) may not be where the stream was (elevation = 0)
    # need to find that point, and adjust the alpha values accordingly 

    # alternatively could use the streamline and find the nearest point
    temp = xs_points.copy()

    offsets = {}
    for ind in temp['cross_section_id'].unique():
        df = temp.loc[temp['cross_section_id'] == ind]

        if df.loc[df['alpha'] == 0]['elevation'].iloc[0] != min(df['elevation']):
            min_ind = df['elevation'].idxmin(skipna=True)        
            temp.loc[temp['cross_section_id'] == ind, 'alpha'] = (df['alpha'] - df['alpha'][min_ind])

    return temp

def find_break_point(df, peak_threshold=0.002, slope_threshold=.4):
    """ df is a half profile dataframe """
    inds, _ =  signal.find_peaks(-df['curvature'], height=peak_threshold)
    for i,ind in enumerate(inds):
        if i == len(inds) - 1:
           slope = (df['elevation'].iloc[-1] - df['elevation'].iloc[ind]) / (df['alpha'].iloc[-1] - df['alpha'].iloc[ind])
        else:
            slope = (df['elevation'].iloc[inds[i+1]] - df['elevation'].iloc[ind]) / (df['alpha'].iloc[inds[i+1]] - df['alpha'].iloc[ind])

        if slope >= slope_threshold:
            return df['point_id'].iloc[ind]
    return 

def determine_threshold(break_points):
    threshold = break_points['hand'].quantile(.6)
    return threshold

def rioxarray_sample_points(raster, points, method='nearest'):
    """ Sample points from raster using rioxarray """
    xs = xr.DataArray(points.geometry.x.values, dims='z')
    ys = xr.DataArray(points.geometry.y.values, dims='z')
    values = raster.sel(x=xs, y=ys, method=method).values
    return values

def add_attributes_to_xs(wbt, points, dem, stream, flow_dir):
    # save rasters to temp files 
    dem_file = os.path.join(wbt.work_dir, 'temp_dem.tif')
    stream_file = os.path.join(wbt.work_dir, 'temp_stream.tif')
    flow_dir_file = os.path.join(wbt.work_dir, 'temp_flow_dir.tif')

    dem.rio.to_raster(dem_file)
    stream.rio.to_raster(stream_file)
    flow_dir.rio.to_raster(flow_dir_file)

    rasters = compute_terrain_rasters(wbt, dem_file, stream_file, flow_dir_file, sigma=1.5)
    dem = rioxarray.open_rasterio(dem_file, masked=True).squeeze()
    slope = rioxarray.open_rasterio(rasters['slope'], masked=True).squeeze()
    curvature = rioxarray.open_rasterio(rasters['curvature'], masked=True).squeeze()
    hand = rioxarray.open_rasterio(rasters['hand'], masked=True).squeeze()

    points['elevation'] = rioxarray_sample_points(dem, points)
    points['slope'] = rioxarray_sample_points(slope, points)
    points['curvature'] = rioxarray_sample_points(curvature, points)
    points['hand'] = rioxarray_sample_points(hand, points)
    points['point_id'] = points.index

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
    
def points_to_polygon(break_points, hand_raster, slope_raster, 
                      slope_threshold=.4, buffer=1):
    threshold = determine_threshold(break_points)

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

def delineate_valley(wbt, dem_raster, flowline, stream_raster, hillslope_raster, flow_dir_raster):
    flowline = flowline.simplify(10)
    points = get_cross_section_points(flowline, xs_spacing=100, xs_width=200, xs_point_spacing=10)

    points = add_attributes_to_xs(wbt, points, dem_raster, stream_raster, flow_dir_raster)

    break_points = []
    for xs_id in points['cross_section_id'].unique():
        df = points.loc[points['cross_section_id'] == xs_id]
        df = rezero_xsection(df)

        # positive break point
        df_pos = df.loc[df['alpha'] >= 0]
        pos_break_point = find_break_point(df_pos)

        # negative break point
        df_neg = df.loc[df['alpha'] <= 0]
        neg_break_point = find_break_point(df_neg)

        break_points.append(pos_break_point)
        break_points.append(neg_break_point)

    break_points = pd.concat(break_points)
    valley_polygon = points_to_polygon(points, hillslope)

    return {'geometry': valley_polygon, 'hand_threshold': threshold}
