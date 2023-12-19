"""
This module takes a dem and creates a subbasin shapefile.

# input:
    dem: dem raster
    threshold: threshold for stream delineation

# output:
    subbasins: subbasins raster
    filled_dem: dem with depressions filled
    d8_pntr: d8 pointer raster
    streams: streams raster

"""
import os

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
            'subbasins': os.path.join(wbt.work_dir, "subbasins.tif"),
            'streams': streams_raster_file}

