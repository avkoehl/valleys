# download flowlines from the National Hydrography Dataset
# use NHDPlus V2 (medium resolution)

# parameters:
#   region: region of interest
#   output: output file name

from pynhd import NHD
import geopandas as gpd

region = '1806'
output = f'../data/{region}_flowlines.geojson'

regions = gpd.read_file('../data/subregions.geojson')
region = regions[regions['huc4'] == region].geometry.iloc[0].bounds # sometimes the region is a multipolygon, so we need to get the bounds

# confirm that filtering before is same as filtering after,
# wouldn't want to only get flowlines that have values for all fields for example
fields = [
    'StreamOrde',
    'StartFlag',
    'TerminalFl',
    'AreaSqKM',
    'LENGTHKM',
    'ArbolateSu',
    'SLOPE',
    'MAXELEVRAW',
    'MINELEVRAW',
    'QA_MA',
    'VA_MA',
    'QA_01',
    'VA_01']

# get flowlines
nhd = NHD('flowline_mr', outfields = fields)
flowlines = nhd.bygeom(region)

# save flowlines
flowlines.to_file(output, driver='GeoJSON')
