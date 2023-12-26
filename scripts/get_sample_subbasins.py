import py3dep
import geopandas as gpd
from pygeohydro import WBD
import rioxarray

odir =  "../data/sample_subbasins/"
catchments_dir = "../data/sample_subbasins/catchment/"
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

# ----- save into subbasins_dir
# create subbasin raster

# ----- save into catchments_dir
# for each catchment save the data
