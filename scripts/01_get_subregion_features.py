from pygeohydro import WBD

path = "../data/subregions.geojson"

wbd = WBD("huc4")
regions = [str(val) for val in range(1801, 1811)]
print('get regions')
regions = wbd.byids("huc4", regions)

print('write to file')
regions.to_file(path, driver="GeoJSON")
