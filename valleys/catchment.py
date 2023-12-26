"""
This module has functions for working with individual catchments. Catchments
come from a subbasin delineated dem.

# input:
    # run code in subbasins.py first

    dem: dem raster file
    streams: stream raster file
    subbasins: basins raster file
    d8_pntr: d8 pointer raster file
"""
import os

import rioxarray
import geopandas as gpd

def wbt_vectorize_stream(wbt, stream_raster_file, d8_pntr_file, ofile):
    wbt.raster_streams_to_vector(stream_raster_file, d8_pntr_file, ofile)
    return os.path.join(wbt.work_dir, ofile)

def wbt_hand(wbt, dem_raster_file, stream_raster_file, ofile):
    wbt.elevation_above_stream(dem_raster_file, stream_raster_file, ofile)
    return os.path.join(wbt.work_dir, ofile)

def subset_raster(raster_file, subbasin_file, basin_id, output_file):
    raster = rioxarray.open_rasterio(raster_file)
    subbasins = rioxarray.open_rasterio(subbasin_file)
    raster = raster.where(subbasins == basin_id)
    raster = raster.dropna(dim="x", how="all")
    raster = raster.dropna(dim="y", how="all")
    raster.rio.to_raster(output_file)
    return output_file

def wbt_profile_curvature(wbt, dem_file, ofile, smooth=False, sigma=3):
    # sigma is for gaussian filter
    if smooth:
        temp_file = os.path.join(wbt.work_dir, "temp.tif")
        wbt.gaussian_filter(dem_file, temp_file, sigma)
        wbt.profile_curvature(temp_file, ofile)
        os.remove(temp_file)
        return ofile

    wbt.profile_curvature(dem_file, ofile)
    return ofile

def get_catchment_data(wbt, dem_file, subbasin_file, streams_file, d8_pntr_file, basin_id):

    cdf = subset_raster(dem_file, subbasin_file, basin_id, os.path.join(wbt.work_dir, f"{basin_id}_dem.tif"))
    csf = subset_raster(streams_file, subbasin_file, basin_id, os.path.join(wbt.work_dir, f"{basin_id}_stream.tif"))
    cpf = subset_raster(d8_pntr_file, subbasin_file, basin_id, os.path.join(wbt.work_dir, f"{basin_id}_d8_pntr.tif"))

    streamline = wbt_vectorize_stream(wbt, csf, cpf, f"{basin_id}_streamline.shp")
    hand = wbt_hand(wbt, cdf, csf, f"{basin_id}_hand.tif")
    curvature = wbt_profile_curvature(wbt, cdf, f"{basin_id}_profile_curvature.tif")

    return {"hand": hand, "streamline": "profile_curvature": curvature, 
            streamline, "basin_id": basin_id, "dem": cdf, "streams": csf, 
            "d8_pntr": cpf}
