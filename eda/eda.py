import geopandas as gpd
import seaborn as sns
import matplotlib.pyplot as plt

flowlines = gpd.read_file('../data/1806_flowlines.geojson')

cat = ['StreamOrde', 'StartFlag', 'TerminalFl']
num = ['AreaSqKM', 'LENGTHKM', 'ArbolateSu', 'SLOPE', 'MAXELEVRAW', 
       'MINELEVRAW', 'QA_MA', 'VA_MA', 'QA_01', 'VA_01']

# maxelev and minelev are in cm
# qa_ma in cfs 
# va_ma in fps 
# qa_01 in cfs
# va_01 in fps


# need to clean the data
# make sure the data types are correct
# make sure the data is in the correct units
# round to 2 decimal places
# figure out values to map to NaN (e.g -9999.00, -9988.00, 0.00)

