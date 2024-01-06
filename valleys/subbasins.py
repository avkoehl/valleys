"""
This module takes a dem and splits it into subbasins. It also extracts the
hillshlopes and streams for each subbasin.

# input:
    dem: dem raster
    <optional> streams: streams raster
    <optional> threshold: threshold for streams raster (needed if streams not provided)

# output:
    subbasins: raster with subbasin id for each pixel
    hillslopes: raster with hillslope id for each pixel, 0 if pixel is a stream
"""
import os

import rioxarray
import pandas as pd

from valleys.utils import chomp_raster

def extract_streams(wbt, dem_raster_file, threshold):

    # wbt.breach_depressions_least_cost(dem_raster_file, "filled_dem.tif", 
    #                                  dist=10, max_cost=None, flat_increment=None,
    #                                  fill=False, min_dist=True)

    wbt.fill_depressions(dem_raster_file, "filled_dem.tif")
    wbt.d8_pointer("filled_dem.tif", "d8_pntr.tif")
    wbt.d8_flow_accumulation("d8_pntr.tif", "flow_acc.tif", pntr=True)
    wbt.extract_streams("flow_acc.tif", "streams.tif", threshold=threshold)

    # cleanup intermediate files
    os.remove(os.path.join(wbt.work_dir, "filled_dem.tif"))
    os.remove(os.path.join(wbt.work_dir, "d8_pntr.tif"))
    os.remove(os.path.join(wbt.work_dir, "flow_acc.tif"))

    return os.path.join(wbt.work_dir, "streams.tif")

def get_subbasins(wbt, dem_raster_file, streams_raster_file):
    #wbt.breach_depressions_least_cost(dem_raster_file, "filled_dem.tif", 
    #                                  dist=10, max_cost=None, flat_increment=None,
    #                                  fill=False, min_dist=True)
    wbt.fill_depressions(dem_raster_file, "filled_dem.tif")
    wbt.d8_pointer("filled_dem.tif", "d8_pntr.tif")
    wbt.subbasins("d8_pntr.tif", streams_raster_file, "subbasins.tif")
    wbt.hillslopes("d8_pntr.tif", streams_raster_file, "hillslopes.tif")

    return {'dem': dem_raster_file,
            'filled_dem': os.path.join(wbt.work_dir, "filled_dem.tif"),
            'd8_pntr': os.path.join(wbt.work_dir, "d8_pntr.tif"),
            'subbasins': os.path.join(wbt.work_dir, "subbasins.tif"),
            'hillslopes': os.path.join(wbt.work_dir, "hillslopes.tif"),
            'streams': streams_raster_file}


def split_into_subbasins(wbt, dem_raster_file, streams_raster_file, output_dir):
    subbasin_results = get_subbasins(wbt, dem_raster_file, streams_raster_file)

    subbasins = rioxarray.open_rasterio(subbasin_results['subbasins']).squeeze()
    subbasins = subbasins.where(subbasins != -32768)

    hillslopes = rioxarray.open_rasterio(subbasin_results['hillslopes']).squeeze()
    hillslopes = hillslopes.where(hillslopes != -32768)

    dem = rioxarray.open_rasterio(dem_raster_file).squeeze()
    dem = dem.where(dem != -32768)

    for subbasin_id in pd.Series(subbasins.values.flatten()).dropna().unique():
        hillslope_ofile = os.path.join(output_dir,
                                       "{}_hillslopes.tif".format(subbasin_id))
        # crop hillslope raster to subbasin and save
        hs = hillslopes.where(subbasins == subbasin_id)
        hs = chomp_raster(hs)
        hs.rio.to_raster(hillslope_ofile)

        # crop dem raster to subbasin and save
        dem_ofile = os.path.join(output_dir, "{}_dem.tif".format(subbasin_id))
        d = dem.where(subbasins == subbasin_id)
        d = chomp_raster(d)
        d.rio.to_raster(dem_ofile)
