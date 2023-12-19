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

def quick_plot(raster_file):
    raster = rioxarray.open_rasterio(raster_file)
    masked = raster.where(raster != -32768)
    masked.plot()
