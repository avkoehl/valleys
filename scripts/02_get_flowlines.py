# download flowlines from the National Hydrography Dataset
# use NHDPlus V2 (medium resolution)

import glob

import geopandas as gpd
import pandas as pd
from pynhd import NHD

output = f'../data/flowlines_all.geojson'
subregions = gpd.read_file('../data/subregions.geojson')

fields = [
    'COMID',
    'StreamOrde',
    'StartFlag',
    'TerminalFl',
    'AreaSqKM',
    'LENGTHKM',
    'ArbolateSu',
    'SLOPE',
    'MAXELEVSMO',
    'MINELEVSMO',
    'QA_MA',
    'VA_MA',
    'QA_01',
    'VA_01']

nhd = NHD('flowline_mr', outfields = fields)

flowlines_by_region = []
for i,row in subregions[['huc4', 'geometry']].iterrows():
    print(f"getting flowlines for {row['huc4']}")
    aoi = row['geometry']

    if aoi.geom_type == 'MultiPolygon':
        aoi = aoi.bounds
    
    flowlines = nhd.bygeom(aoi)
    flowlines['huc4'] = row['huc4']
    flowlines_by_region.append(flowlines)

flowlines = gpd.GeoDataFrame(pd.concat(
    flowlines_by_region, ignore_index=True), crs=flowlines_by_region[0].crs)
flowlines.to_file(output, driver='GeoJSON')


