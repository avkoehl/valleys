"""
code for finding breakpoints on cross sections and then delineating valley floor from those breakpoints

- dataset: ['smoothed_dem', 'slope', 'curvature', 'streams', 'hillslopes', 'flow_dir', 'hand']
- flowline
-------
- cross_sections_df
- break_points_df
- hand_threshold
- valley_floor_polygon
- valley_floor_raster


sample_cross_section_points
find_breakpoints
determine_hand_threshold
delineate_valley_floor
valley_floor_full_workflow

"""
import os

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import scipy
from shapely.geometry import shape
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
import xarray as xr

from valleys.cross_section import get_cross_section_points
from valleys.breakpoints import find_xs_break_points

class Subbasin:
    def __init__(self, dataset, flowline, subbasin_id):
        self.subbasin_id = subbasin_id
        self.flowline = flowline
        self.dataset = dataset

        required_datasets = ['elevation', 'slope', 'curvature', 'strm_val', 'hillslopes', 'flow_dir', 'hand']
        for band in required_datasets:
            if band not in self.dataset:
                raise ValueError(f"Required band '{band}' not found in the input dataset")

        self.cross_sections_df = None
        self.break_points_df = None
        self.hand_threshold = None
        self.valley_floor_polygon = None
        self.valley_floor_raster = None

    def _preprocess_flowline(self):
        self.flowline = self.flowline.simplify(20)
    
    def sample_cross_section_points(self):
        self._preprocess_flowline()

        points = get_cross_section_points(self.flowline, xs_spacing=20, xs_width=500, xs_point_spacing=10)
        points['point_id'] = np.arange(len(points))

        for data_layer in self.dataset.data_vars:
            points[data_layer] = _rioxarray_sample_points(self.dataset[data_layer], points)

        points = points.loc[~points['elevation'].isna()]
        points = points.loc[~points['slope'].isna()]
        self.cross_sections_df = points

    def find_breakpoints(self):
        break_points_list = []
        for xs in self.cross_sections_df['cross_section_id'].unique():
            xs_points = self.cross_sections_df.loc[self.cross_sections_df['cross_section_id'] == xs]
            break_points = find_xs_break_points(xs_points)
            break_points = (xs, *break_points)
            break_points_list.append(break_points)
        break_points_df = pd.DataFrame(break_points_list, columns=['cross_section_id', 'pos', 'neg', 'peak_ids'])
        combined = break_points_df['pos'].dropna().to_list() + break_points_df['neg'].dropna().to_list()
        self.break_points_df = self.cross_sections_df.loc[self.cross_sections_df['point_id'].isin(combined)]

    def determine_hand_threshold(self):
        hand_values = self.break_points_df['hand']

        # remove outliers
        hand_values = hand_values[hand_values <  hand_values.quantile(.95)]
        hand_values = hand_values[hand_values <  50]

        # set threshold
        self.hand_threshold = hand_values.quantile(.7)

    def delineate_valley_floor(self):
        hand = self.dataset['hand']
        threshold = self.hand_threshold

        values = _apply_threshold_and_fill_holes(hand, threshold)
        values = _combine_with_slope_threshold(values, self.dataset['slope'], 30)
        
        # polygonize
        polygons = _polygonize(values)
        polygons = gpd.GeoDataFrame(geometry=polygons, crs=3310)
        polygons['geometry'] = polygons['geometry'].apply(_close_holes)

        # convert to multipolygon or single polygon
        if len(polygons) > 1:
            polygon = MultiPolygon(polygons['geometry'].values)
            self.valley_floor_polygon = polygon

        if len(polygons) == 1:
            self.valley_floor_polygon = polygons['geometry'].iloc[0]
            self.valley_floor_polygon = polygons['geometry'].iloc[0]
            self.valley_floor_raster = values

        return

    def valley_floor_by_breakpoints_full_workflow(self):
        self.sample_cross_section_points()
        self.find_breakpoints()
        self.determine_hand_threshold()
        self.delineate_valley_floor()
        pass

    def valley_floor_slope_threshold_workflow(self):
        # just use very high hand threshold and trim by slope
        # TODO
        pass

    def plot_hand_mean_slope_relation(self, odir):
        # TODO
        # try this on whole watershed
        pass

    def plot_hand_cdf(self, odir):
        # is self.hand_threshold add that as a vertical line
        # else plot as normal
        # try this on whole watershed
        # TODO
        pass

    def plot_breakpoints(self, odir):
        # TODO 
        # create matplotlib figure with cross section points and breakpoints and peaks and 'channel' for each cross section
        # and dump into folder (odir)
        pass

def _combine_with_slope_threshold(binary_raster, slope_raster, slope_threshold):
    values = binary_raster.where((slope_raster < slope_threshold) & (binary_raster == 1))
    values.data = values.data.astype(np.uint8)
    values.data = scipy.ndimage.binary_fill_holes(values.data)
    values = values.where(values != 0)
    return values


def _apply_threshold_and_fill_holes(raster, threshold):
    values = raster.where(raster < threshold)
    values = values.where(np.isnan(values), 1)
    values = values.where(~np.isnan(values), 0)
    values = (values > 0).astype(int)

        # fill
    values.data = scipy.ndimage.binary_fill_holes(values.data)
    values = values.where(values != 0)
    return values


def _rioxarray_sample_points(raster, points, method='nearest'):
    xs = xr.DataArray(points.geometry.x.values, dims='z')
    ys = xr.DataArray(points.geometry.y.values, dims='z')
    values = raster.sel(x=xs, y=ys, method=method).values
    return values

def _close_holes(poly):
    if len(poly.interiors):
        return Polygon(list(poly.exterior.coords))
    return poly

def _polygonize(raster):
    # binary raster 1,0
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
