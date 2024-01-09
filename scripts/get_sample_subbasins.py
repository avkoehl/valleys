import glob
import os
import shutil

import py3dep
import numpy as np
import pandas as pd
import geopandas as gpd
from pygeohydro import WBD
import rioxarray

from valleys.subbasins import extract_streams
from valleys.subbasins import split_into_subbasins
from valleys.utils import setup_wbt

odir = "../data/huc12s/"

if os.path.exists(odir):
    shutil.rmtree(odir)
os.mkdir(odir)

wbd = WBD("huc12")
wbt = setup_wbt(whitebox_dir="~/opt/WBT/", working_dir="../data/wb_outputs/")

# get a list of huc12s that cover a variety of terrains
huc12s = [
        '180701020604',
        '180101100904',
        '180400060403']

for huc12 in huc12s:
    print("processing huc12: ", huc12)
    huc12_dir = f"{odir}{huc12}/"
    if os.path.exists(huc12_dir):
        shutil.rmtree(huc12_dir)
    os.mkdir(huc12_dir)
    os.mkdir(f"{huc12_dir}/subbasins/")

    boundary = wbd.byids("huc12", huc12).iloc[0]
    dem = py3dep.get_map("DEM", boundary.geometry, resolution=10,
                         geo_crs="epsg:4326", crs="4326")
    dem = dem.rio.reproject(3310)
    dem.rio.to_raster(f"{huc12_dir}{huc12}_dem_10m.tif")

    stream_file = extract_streams(wbt, 
                        os.path.abspath(f"{huc12_dir}/{huc12}_dem_10m.tif"), 
                                    threshold=20000)

    split_into_subbasins(wbt, os.path.abspath(f"{huc12_dir}{huc12}_dem_10m.tif"),
                         stream_file, os.path.abspath(f"{huc12_dir}/subbasins/"))
