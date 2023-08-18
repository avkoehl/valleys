import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

flowlines = gpd.read_file('../data/flowlines_all.geojson')

cat = ['StreamOrde', 'StartFlag', 'TerminalFl']
num = ['AreaSqKM', 'LENGTHKM', 'ArbolateSu', 'SLOPE', 'MAXELEVRAW', 
       'MINELEVRAW', 'QA_MA', 'VA_MA', 'QA_01', 'VA_01']

# Data Cleaning ----------------------------------------------------------------
# convert all NHD encoded missing values to NaN
flowlines = flowlines.replace(to_replace=[-9999.00, -9998.00], value=np.nan)

for col in ['QA_MA', 'VA_MA', 'QA_01', 'VA_01']:
    flowlines[col] = flowlines[col].replace(to_replace=0.00, value=np.nan)

for col in ['AreaSqKM', 'LENGTHKM', 'ArbolateSu','MAXELEVRAW', 'MINELEVRAW']:
    flowlines[col] = flowlines[col].replace(to_replace=0.00, value=np.nan)

# if slope is greater than 1, set to NaN
flowlines['SLOPE'].where(flowlines['SLOPE'] < 1.00, np.nan, inplace=True)

# if stream order is < 1, set to NaN
flowlines['StreamOrde'].where(flowlines['StreamOrde'] > 0, np.nan, inplace=True)

# convert units
flowlines['MAXELEVSMO'] = flowlines['MAXELEVSMO'] / 100 # cm to m
flowlines['MINELEVSMO'] = flowlines['MINELEVSMO'] / 100 # cm to m

# Data Exploration -------------------------------------------------------------
# get counts of missing values for each column
(flowlines.isna().sum() / len(flowlines) * 100).round()

# count plots of categorical variables X SubRegion
fig, ax = plt.subplots(1,3, figsize=(10, 10))
for i,col in enumerate(flowlines[cat]):
    sns.countplot(x=col, data=flowlines, ax=ax[i])

# histogram of numerical variables

# pairplot of numerical variables

# boxplots of numerical variables X categorical variables

