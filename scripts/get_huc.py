import os

from pygeohydro import WBD
from pynhd import NHD
import py3dep

def get_dem_and_flowlines(hucid, layer):
    wbd = WBD(layer)
    boundary = wbd.byids(layer, hucid)
    boundary_reprojected = boundary.to_crs(3310)

    dem = py3dep.get_map("DEM", boundary.geometry.iloc[0], resolution=10, geo_crs=boundary.crs, crs=4326)
    dem = dem.rio.reproject(3310)

    nhd_mr = NHD("flowline_mr")
    flowlines_mr = nhd_mr.bygeom(boundary.geometry.iloc[0].bounds)
    flowlines_mr = flowlines_mr.to_crs(3310)
    flowlines_mr = flowlines_mr.clip(boundary_reprojected.geometry.iloc[0])

    return dem, flowlines_mr

def get_huc(hucid, layer, odir="../sampledata/"):
    odir = os.path.join(odir, hucid)
    if not os.path.exists(odir):
        os.makedirs(odir)
    dem, flowlines_mr = get_dem_and_flowlines(hucid, layer)
    dem.rio.to_raster(os.path.join(odir, "dem.tif"))
    flowlines_mr.to_file(os.path.join(odir, "flowlines_mr.shp"))

#huc12s = ["180600060101", "180400060301", "180701060201"]
#for huc12 in huc12s:
#    print(huc12)
#    get_huc(huc12, 'huc12')

huc10s = ["1801010701", "1802012804", "1803001003"]
for huc10 in huc10s:
    print(huc10)
    get_huc(huc10, 'huc10', odir='../data/')

