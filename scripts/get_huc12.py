import os

from pygeohydro import WBD
from pynhd import NHD
import py3dep

def get_dem_and_flowlines(huc12):
    wbd = WBD("huc12")
    boundary = wbd.byids("huc12", huc12)
    boundary_reprojected = boundary.to_crs(3310)

    dem = py3dep.get_map("DEM", boundary.geometry.iloc[0], resolution=10, geo_crs=boundary.crs, crs=4326)
    dem = dem.rio.reproject(3310)

    nhd_mr = NHD("flowline_mr")
    flowlines_mr = nhd_mr.bygeom(boundary.geometry.iloc[0].bounds)
    flowlines_mr = flowlines_mr.to_crs(3310)
    flowlines_mr = flowlines_mr.clip(boundary_reprojected.geometry.iloc[0])

    return dem, flowlines_mr

huc12 = "180600060101"
odir = f"../sampledata/{huc12}/"

dem, flowlines_mr = get_dem_and_flowlines(huc12)

if not os.path.exists(odir):
    os.makedirs(odir)

dem.rio.to_raster(os.path.join(odir, "dem.tif"))
flowlines_mr.to_file(os.path.join(odir, "flowlines_mr.shp"))
