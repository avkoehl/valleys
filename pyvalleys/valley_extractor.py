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
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.geometry import MultiPoint
from shapely.geometry import MultiPolygon
from shapelysmooth import taubin_smooth
from shapelysmooth import chaikin_smooth
import xarray as xr

from pyvalleys.cross_section import get_cross_section_points
from pyvalleys.cross_section import get_cross_section_points_from_points
from pyvalleys.cross_section import get_cross_section_lines_from_points
from pyvalleys.breakpoints import find_xs_break_points
from pyvalleys.breakpoints import find_xs_break_points_alternate
from pyvalleys.gis import rioxarray_sample_points, close_holes, polygonize_feature, get_length_and_width, get_points_on_linestring

class ValleyExtractor:
    def __init__(self, dataset, flowline, subbasin_id):
        self.subbasin_id = subbasin_id
        self.flowline = flowline
        self.flowline_raw = flowline
        self.dataset = dataset

        required_datasets = ['elevation', 'slope', 'curvature', 'strm_val', 'hillslopes', 'flow_dir', 'hand']
        for band in required_datasets:
            if band not in self.dataset:
                raise ValueError(f"Required band '{band}' not found in the input dataset")

        self.cross_sections_df = None
        self.cross_sections_lines = None
        self.break_points_df = None
        self.hand_threshold = None
        self.valley_floor_polygon = None
        self.valley_floor_raster = None
        self.polygon = None

    def get_region_polygon(self):
        # hand most likely to have issues
        elevation = self.dataset['elevation']
        values = elevation
        values = values.where(np.isnan(values), 1)
        values = values.astype(np.uint8)

        polygons = polygonize_feature(values, 1)
        if len(polygons) == 0:
            raise RuntimeError("something wrong polygonizing the hand raster")
        if len(polygons) == 1:
            polygon = polygons[0]
        if len(polygons) > 1:
            series = gpd.GeoSeries(polygons).buffer(1)

            if len(series) > 1:
                # remove tiny polygons
                series = series.loc[series.area > 1000]
                if len(series) > 1:
                    polygon = MultiPolygon(series.values)
                else:
                    polygon = series.iloc[0]
            else:
                polygon = polygon.iloc[0]

        self.polygon = polygon
    
    def _get_max_possible_width(self):
        if self.polygon is None:
            self.get_region_polygon()
        return max(get_length_and_width(self.polygon))

    def _preprocess_flowline(self, tolerance=None, smooth=True):
        if tolerance:
            self.flowline = self.flowline.simplify(tolerance)
        else:
            self.flowline = self.flowline

        if smooth:
            self.flowline = chaikin_smooth(taubin_smooth(self.flowline))

    def get_cross_section_data(self, tolerance, smooth, xs_spacing, xs_width, xs_point_spacing):
        # get cross sections:
        #      1. cross_section_id, point_id, point, value1, value2 ...
        #      3. cross_section_id, linestring (clipped to subbasin)
        self._preprocess_flowline(tolerance=tolerance, smooth=smooth)

        points = get_points_on_linestring(self.flowline, xs_spacing)

        # 1.
        xs_points = get_cross_section_points_from_points(self.flowline, points, xs_width, xs_point_spacing)

        # 2.
        max_width = self._get_max_possible_width()
        width = int(max_width + 1) #this will completely overshoot but that's probably okay
        lines = get_cross_section_lines_from_points(self.flowline, points, width)
        lines = lines.clip(self.polygon)
        lines = lines.explode(index_parts=False)
        lines = lines.loc[lines.intersects(self.flowline.buffer(1))]

        # TODO: sometime flowline is outside of polygon slightly
        self.cross_sections_lines = lines


        xs_points['point_id'] = np.arange(len(xs_points))
        for data_layer in self.dataset.data_vars:
            xs_points[data_layer] = rioxarray_sample_points(self.dataset[data_layer], xs_points)
        xs_points = xs_points.loc[~xs_points['elevation'].isna()]
        xs_points = xs_points.loc[~xs_points['slope'].isna()]
        self.cross_sections_df = xs_points
        return
    
    def sample_cross_section_points(self, tolerance, xs_spacing, xs_width, xs_point_spacing):
        self._preprocess_flowline(tolerance)

        points = get_cross_section_points(self.flowline, xs_spacing=xs_spacing, xs_width=xs_width, xs_point_spacing=xs_point_spacing)
        points['point_id'] = np.arange(len(points))

        for data_layer in self.dataset.data_vars:
            points[data_layer] = rioxarray_sample_points(self.dataset[data_layer], points)

        points = points.loc[~points['elevation'].isna()]
        points = points.loc[~points['slope'].isna()]
        self.cross_sections_df = points

    def find_breakpoints(self, peak_threshold=0.002, bp_slope_threshold=20, method=None):
        break_points_list = []
        for xs in self.cross_sections_df['cross_section_id'].unique():
            xs_points = self.cross_sections_df.loc[self.cross_sections_df['cross_section_id'] == xs]

            if method == 'alternate':
                # this use curvature of the cross section profile not the profile of curvatures
                break_points = find_xs_break_points_alternate(xs_points, peak_threshold=peak_threshold, slope_threshold=bp_slope_threshold)
            else:
                break_points = find_xs_break_points(xs_points, peak_threshold=peak_threshold, slope_threshold=bp_slope_threshold)

            break_points = (xs, *break_points)
            break_points_list.append(break_points)
        break_points_df = pd.DataFrame(break_points_list, columns=['cross_section_id', 'pos', 'neg', 'peak_ids'])
        combined = break_points_df['pos'].dropna().to_list() + break_points_df['neg'].dropna().to_list()
        self.break_points_df = self.cross_sections_df.loc[self.cross_sections_df['point_id'].isin(combined)]

    def determine_hand_threshold(self, quantile=0.7, buffer=0):
        hand_values = self.break_points_df['hand']

        # remove outliers
        hand_values = hand_values[hand_values <  hand_values.quantile(.95)]
        #hand_values = hand_values[hand_values <  50]

        # set threshold
        self.hand_threshold = hand_values.quantile(quantile)

        if buffer > 0:
            self.hand_threshold = self.hand_threshold + buffer

    def delineate_valley_floor(self, slope_threshold=None):
        values = self._apply_thresholds_and_fill_holes(self.hand_threshold, slope_threshold)
        self.valley_floor_raster = values
        
        # polygonize
        # binarize to 0 and 1
        values = (values > 0).astype(np.uint8)
        polygons = polygonize_feature(values, 1) # get polygons where raster value is 1
        polygons = [close_holes(p) for p in polygons]
        polygons = gpd.GeoDataFrame(geometry=polygons, crs=3310)
        # remove polygons that dont intersect the flowline
        polygons = polygons.loc[polygons.intersects(self.flowline_raw)]

        # convert to multipolygon or single polygon
        if len(polygons) > 1:
            polygon = MultiPolygon(polygons['geometry'].values)
            polygon = polygon.buffer(1)  # cleans up multipolygon that is just an artefact of the raster to vector
            polygon = close_holes(polygon)

            self.valley_floor_polygon = polygon.iloc[0]

        if len(polygons) == 1:
            self.valley_floor_polygon = polygons['geometry'].iloc[0]

        return

    def _apply_thresholds_and_fill_holes(self, hand_threshold, slope_threshold=25):
        values = self.dataset['hand'].where(self.dataset['hand'] <= hand_threshold)
        # binarize
        values = values.where(np.isnan(values), 1)
        values = values.where(~np.isnan(values), 0)
        values = (values > 0).astype(np.uint8)
        # fill holes
        values.data = scipy.ndimage.binary_fill_holes(values.data)
    
        # Optional apply slope threshold
        if slope_threshold:
            slope = self.dataset['slope']
            values = values.where(slope <= slope_threshold)
            values = (values > 0).astype(np.uint8)
            # fill holes again
            values.data = scipy.ndimage.binary_fill_holes(values.data)
    
        # burnin streams
        values = values.where(self.dataset['strm_val'] != self.subbasin_id, 1)
    
        # assign subbasin id to valley floor
        values = values.where(values != 1, self.subbasin_id)
        return values

    def run(self, smooth=True, tolerance=20, xs_spacing=20,
                                                  xs_width=500, xs_point_spacing=10,
                                                  quantile=0.7, buffer=0, 
                                                  slope_threshold=None, peak_threshold=0.002, bp_slope_threshold=20):

        NOT_ENOUGH_WALLS = False
        self.get_cross_section_data(tolerance=tolerance, smooth=smooth, xs_spacing=xs_spacing, xs_width=xs_width, xs_point_spacing=xs_point_spacing)
        self.find_breakpoints(peak_threshold=peak_threshold, bp_slope_threshold=bp_slope_threshold)

        # if not enough breakpoints: rerun with xs_width == max width of subbasin
        if len(self.break_points_df) < 5:
            new_xs_width = int(self._get_max_possible_width() + 1)
            self.get_cross_section_data(tolerance=tolerance, smooth=smooth, xs_spacing=xs_spacing, xs_width=new_xs_width, xs_point_spacing=xs_point_spacing)
            self.find_breakpoints(peak_threshold=peak_threshold, bp_slope_threshold=bp_slope_threshold)

            # if still not enough
            if len(self.break_points_df) < 5:
                # self.hand_threshold = self.dataset['hand'].max().values.item()
                points = []
                for line in self.cross_sections_lines['geometry']:
                    start = Point(line.coords[0])
                    end = Point(line.coords[-1])
                    points.append(start)
                    points.append(end)

                points = MultiPoint(points)
                polygon = points.convex_hull
                polygon = gpd.GeoSeries(polygon).clip(self.polygon)
                NOT_ENOUGH_WALLS = True

            else:
                self.determine_hand_threshold(quantile, buffer)
        else:
            self.determine_hand_threshold(quantile, buffer)
            
        if not NOT_ENOUGH_WALLS:
            self.delineate_valley_floor(slope_threshold=slope_threshold)
        else:
            print("no walls found by cross sections")
            self.hand_threshold = None
            self.valley_floor_polygon = polygon
            # alteranatively set threshold based on quantile of HAND values of boundary points?
            # still need to apply slope threshold and clip better
        return
