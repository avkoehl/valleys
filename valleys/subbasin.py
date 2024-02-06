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
import numpy as np
import xarray as xr

from valleys.cross_section import get_cross_section_points
from valleys.breakpoints import find_xs_break_point

class Subbasin:
    def __init__(self, dataset, flowline, subbasin_id):
        self.subbasin_id = subbasin_id
        self.dataset = dataset
        self.flowline = flowline

        required_datasets = ['elevation', 'slope', 'curvature', 'strm_val', 'hillslopes', 'flow_dir', 'hand']
        for band in required_datasets:
            if band not in self.dataset:
                raise ValueError(f"Required band '{band}' not found in the input dataset")

        self.cross_sections_df = None
        self.break_points_df = None
        self.hand_threshold = None
        self.valley_floor_polygon = None
        self.valley_floor_raster = None
    
    def sample_cross_section_points(self):
        points = get_cross_section_points(self.flowline, xs_spacing=20, xs_width=500, xs_point_spacing=10)
        points['point_id'] = np.arange(len(points))

        for data_layer in self.dataset.data_vars:
            points[data_layer] = _rioxarray_sample_points(self.dataset[data_layer], points)

        points = points.loc[~points['elevation'].isna()]
        points = points.loc[~points['slope'].isna()]
        self.cross_sections_df = points

    def find_breakpoints(self):
        break_points = []
        for xs in self.cross_sections_df['cross_section_id'].unique():
            xs_points = self.cross_sections_df.loc[self.cross_sections_df['cross_section_id'] == xs]
            break_point = find_xs_break_point(xs_points)
            break_points.append(break_point)
        break_points_df = pd.DataFrame(break_points, columns=['cross_section_id', 'pos', 'neg', 'peak_ids'])
        combined = break_points_df['pos'].dropna().to_list() + break_points_df['neg'].dropna().to_list()
        self.break_points_df = self.cross_sections_df.loc[self.cross_sections_df['point_id'].isin(combined)]

    def determine_hand_threshold(self):
        self.hand_threshold = self.break_points_df['hand'].quantile(.8)

    def delineate_valley_floor(self):
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
            self.valley_floor_polygon = polygon
        self.valley_floor_polygon = polygons['geometry'].iloc[0]
        self.valley_floor_raster = values
        pass

    def valley_floor_full_workflow(self):
        self.sample_cross_section_points()
        self.find_breakpoints()
        self.determine_hand_threshold()
        self.delineate_valley_floor()
        pass

def _rioxarray_sample_points(data, points):
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
