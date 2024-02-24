import rasterio
from rasterio import features
import xarray as xr
from shapely.geometry import shape, Polygon

def rioxarray_sample_points(raster, points, method='nearest'):
    xs = xr.DataArray(points.geometry.x.values, dims='z')
    ys = xr.DataArray(points.geometry.y.values, dims='z')
    values = raster.sel(x=xs, y=ys, method=method).values
    return values

def close_holes(poly):
    if len(poly.interiors):
        return Polygon(list(poly.exterior.coords))
    return poly

def polygonize_feature(raster, feature_value=1):
    # get polygons that correspond to the regions where raster == feature_value
    # wrapper for rasterio.features.shapes
    # raster is a numeric array and feature_value is the value of the feature to polygonize

    # create a mask where the raster is equal to the feature value
    # this will be used to mask the raster when calling rasterio.features.shapes
    # because we only want to have to operate on the pixels that are equal to the feature value
    mask = raster == feature_value
        
    polygons = []
    for geom, value in features.shapes(raster, mask=mask, transform=raster.rio.transform()):
        if value == feature_value:  #
            # load the geometry as a shapely Polygon and append to the list
           polygons.append(shape(geom))
    return polygons

def chomp_raster(raster):
    raster = raster.dropna(dim='x', how='all')
    raster = raster.dropna(dim='y', how='all')
    return raster
