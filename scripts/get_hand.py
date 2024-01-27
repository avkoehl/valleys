"""
https://www.hydroshare.org/resource/73aaa3efcda2465ba6227f535400f36b/
"""


import os
import sys

import geopandas as gpd
import pandas as pd
import pynhd
from pygeohydro import WBD
import py3dep
import pygeohydro
import rioxarray as rxr 
import requests
import zipfile



huc12code = 180600060101
huc6code = 180600

odir = "../data/hand/"

if os.path.exists(odir):
    shutil.rmtree(odir)
os.makedirs(odir)

wbd = WBD("huc12")
huc12_boundary =  wbd.byids("huc12", huc12code)

#
# download the huc6 data
# unzip
# clip to huc12 boundary

# https://cfim.ornl.gov/data/HAND/20200601/{huc6code}.zip
url = f"https://cfim.ornl.gov/data/HAND/20200601/{huc6code}.zip"
r = requests.get(url)

with open(f"{odir}/{huc6code}.zip", "wb") as f:
    f.write(r.content)

# unzip
with zipfile.ZipFile(f"{odir}/{huc6code}.zip", 'r') as zip_ref:
    zip_ref.extractall(f"{odir}/")


