import glob
import os
import shutil

import py3dep
import numpy as np
import pandas as pd
import geopandas as gpd
from pygeohydro import WBD
import rioxarray

from valleys.catchment import get_catchment_data
from valleys.subbasins import delineate_streams
from valleys.subbasins import delineate_subbasins
from valleys.utils import setup_wbt


subbasins_dir = "../data/sample_subbasins/subbasin/"
huc12_dir = "../data/sample_subbasins/huc_12/"


# get a list of huc12s that cover a variety of terrains
huc12s = [
        '180701020604',
        '180101100904',
        '180400060403']

# ----- save into huc12_dir
# get polygons for each huc12
wbd = WBD("huc12")
boundaries = wbd.byids("huc12", huc12s)

# for each polygon get the dem
for index,boundary in boundaries.iterrows():
    print(f"Getting DEM for {boundary['huc12']}")
    dem = py3dep.get_map("DEM", boundary.geometry, resolution=10, geo_crs="epsg:4326", crs="4326")
    # convert to 3310 crs
    dem = dem.rio.reproject(3310)
    dem.rio.to_raster(f"{huc12_dir}{boundary['huc12']}_dem_10m.tif")

boundaries.to_file(f"{huc12_dir}boundaries.geojson", driver="GeoJSON")

# setup wbt
wbt = setup_wbt(whitebox_dir="~/opt/WBT/", working_dir="../data/wb_outputs/")

# ----- save into subbasins_dir
# create subbasin raster
# iterate through the dems
records = [] 
for dem_file in glob.glob(f"{huc12_dir}*_dem_10m.tif"):
    print(f"Creating subbasin raster for {dem_file}")

    # delineate streams
    streams_file = delineate_streams(wbt, os.path.abspath(dem_file), threshold=20000)

    # delineate subbasins
    subbasins_info = delineate_subbasins(wbt, os.path.abspath(dem_file), streams_file)
    # keys: filled_dem, d8_pntr, subbasins, streams

    # move files to subbasins_dir for that huc12
    # create dir for that huc12
    huc12 = dem_file.split("/")[-1].split("_")[0]
    if os.path.exists(f"{subbasins_dir}{huc12}"):
        shutil.rmtree(f"{subbasins_dir}{huc12}")
    os.mkdir(f"{subbasins_dir}{huc12}")
    # move files
    for key in subbasins_info.keys():
        os.rename(subbasins_info[key], f"{subbasins_dir}{huc12}/{huc12}_{key}.tif")
        subbasins_info[key] = f"{subbasins_dir}{huc12}/{huc12}_{key}.tif"
        records.append({"huc12": huc12, "key": key, "path": subbasins_info[key]})

df = pd.DataFrame.from_records(records)
df.columns = ["huc12", "raster_name", "path"]
df["path"] = df["path"].apply(lambda x: os.path.abspath(x))

# save df
df.to_csv(f"{subbasins_dir}subbasins.csv", index=False)
