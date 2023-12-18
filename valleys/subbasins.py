"""
This module takes a dem and creates a subbasin shapefile.

# input:
    dem: dem raster
    threshold: threshold for stream delineation

# output:
    subbasins: shapefile of subbasins with column that has the stream


Example usage:
# --- Set up --- #
wbt = setup_wbt("~/opt/WBT/", "../data/whitebox_outputs/")
dem_file = os.path.abspath("../toy_data/dem_3m.tif")
threshold = 30000

# --- Extract Subbasins --- #
streams_raster_file = delineate_streams(wbt, dem_file, threshold)
basins_results = delineate_subbasins(wbt, dem_file, streams_raster_file)

# --- Vectorize Subbasins --- #
subbasins = vectorize_and_merge(wbt, basins_results['subbasins'], 
                                streams_raster_file, basins_results['d8_pntr'])
"""
import os
import shutil

import geopandas as gpd
import rioxarray
import rasterio
from rasterio import features
import shapely
import whitebox

def delineate_streams(wbt, dem_raster_file, threshold):
    wbt.breach_depressions_least_cost(dem_raster_file, "filled_dem.tif", dist=1000)
    wbt.d8_pointer("filled_dem.tif", "d8_pntr.tif")
    wbt.d8_flow_accumulation("d8_pntr.tif", "flow_acc.tif", pntr=True)
    wbt.extract_streams("flow_acc.tif", "streams.tif", threshold=threshold)

    # cleanup intermediate files
    os.remove(os.path.join(wbt.work_dir, "filled_dem.tif"))
    os.remove(os.path.join(wbt.work_dir, "d8_pntr.tif"))
    os.remove(os.path.join(wbt.work_dir, "flow_acc.tif"))

    return os.path.join(wbt.work_dir, "streams.tif")

def delineate_subbasins(wbt, dem_raster_file, streams_raster_file):
    wbt.breach_depressions_least_cost(dem_raster_file, "filled_dem.tif", dist=1000)
    wbt.d8_pointer("filled_dem.tif", "d8_pntr.tif")
    wbt.subbasins("d8_pntr.tif", streams_raster_file, "subbasins.tif")

    return {'filled_dem': os.path.join(wbt.work_dir, "filled_dem.tif"),
            'd8_pntr': os.path.join(wbt.work_dir, "d8_pntr.tif"),
            'subbasins': os.path.join(wbt.work_dir, "subbasins.tif")}

def setup_wbt(whitebox_dir, working_dir):
    wbt = whitebox.WhiteboxTools()
    wbt.set_whitebox_dir(os.path.expanduser(whitebox_dir))
    
    working_directory = os.path.abspath(working_dir)
    if os.path.exists(working_directory):
            shutil.rmtree(working_directory)
            os.mkdir(working_directory)
    wbt.set_working_dir(os.path.abspath(working_directory))
    wbt.set_verbose_mode(False)
    return wbt

def vectorize_subbasins(subbasins_raster_file):
    with rasterio.open(subbasins_raster_file) as src:
        subbasins = src.read(1)
    mask = subbasins != -32768
    shapes = features.shapes(subbasins, mask=mask, transform=src.transform)

    records = []
    for shape,id in shapes:
        coords = shape["coordinates"][0]
        polygon = shapely.Polygon(coords)
        records.append({"geometry": polygon, "subbasin_id": id})
    subbasins = gpd.GeoDataFrame.from_records(records)
    subbasins.crs = src.crs

    subbasins = subbasins[subbasins.area > 100] # remove artifacts
    if subbasins['subbasin_id'].nunique() != len(subbasins):
        print("Warning: subbasins have duplicate ids")

    subbasins.reset_index(inplace=True, drop=True)
    return subbasins

def vectorize_streams(wbt, streams_raster_file, d8_pntr_file):
    wbt.raster_streams_to_vector(streams_raster_file, d8_pntr_file, 
                                 "streams.shp")
    streams = gpd.read_file(os.path.join(wbt.work_dir, "streams.shp"))
    streams.crs = rasterio.open(streams_raster_file).crs
    return streams

def vectorize_and_merge(wbt, subbasins_raster_file, streams_raster_file,
                        d8_pntr_file):
    subbasins = vectorize_subbasins(subbasins_raster_file)
    streams = vectorize_streams(wbt, streams_raster_file, d8_pntr_file)
    # associate subbasin id with stream
    streams['representative_coords'] = streams['geometry'].apply(
            lambda x: x.representative_point())
    streams['subbasin_id'] = streams['representative_coords'].apply(
            lambda x: subbasins['subbasin_id'][subbasins.contains(x)].values[0])
    # confirm that all streams have a unique subbasin id
    assert len(streams) == len(streams['subbasin_id'].unique())
    # confirm that all streams have a subbasin id
    assert len(streams[streams['subbasin_id'].isna()]) == 0
    
    # merge 
    subbasins = subbasins.merge(streams[['subbasin_id', 'geometry']], 
                                on='subbasin_id', how='left')
    subbasins.rename(
            columns={'geometry_x': 'geometry', 'geometry_y': 'stream_geometry'}, 
            inplace=True)
    subbasins = subbasins[['subbasin_id', 'geometry', 'stream_geometry']]
    subbasins = gpd.GeoDataFrame(subbasins, geometry='geometry', 
                                 crs=streams.crs)

    subbasins['stream_geometry'] = subbasins.apply(lambda row:
            row['stream_geometry'].intersection(row['geometry']), axis=1)
    return subbasins

def quick_plot(raster_file):
    raster = rioxarray.open_rasterio(raster_file)
    masked = raster.where(raster != -32768)
    masked.plot()

wbt = setup_wbt("~/opt/WBT/", "../data/whitebox_outputs/")
dem_file = os.path.abspath("../toy_data/dem_3m.tif")
threshold = 30000

# --- Extract Subbasins --- #
streams_raster_file = delineate_streams(wbt, dem_file, threshold)
basins_results = delineate_subbasins(wbt, dem_file, streams_raster_file)
