"""
180101100904     basinid 12.    
180701020604     basin1d 4     
180400060403    basinid 22
bear creek     basinid 10


crop dem to basinid
crop streams to basinid

# confirm that only 1 stream

"""
import os
import shutil

import numpy as np
import pandas as pd
import geopandas as gpd

from valleys.catchment import subset_raster

# make new output folder 
odir = "../toydata/example_valleys/"
if os.path.exists(odir):
    shutil.rmtree(odir)
os.makedirs(odir)

# paths to dems and streams
paths_df = pd.DataFrame(
        {'region': ['180101100904', '180701020604', '180400060403'],
         'basinid': [12, 4, 22]})
paths_df['subbasin'] = paths_df['region'].apply(lambda x: f"../data/sample_subbasins/subbasin/{x}/{x}_subbasins.tif")
paths_df['dem'] = paths_df['region'].apply(lambda x: f"../data/sample_subbasins/subbasin/{x}/{x}_dem.tif")
paths_df['stream_ras'] = paths_df['region'].apply(lambda x: f"../data/sample_subbasins/subbasin/{x}/{x}_streams.tif")

# bear creek
bear_creek = pd.DataFrame({'region': ['bear_creek'], 'basinid': [10], 'subbasin': ['../toy_data/subbasins.tif'], 'dem': ['../toy_data/dem_3m.tif'], 'stream_ras': ['../toy_data/streams.tif']})
paths_df = pd.concat([paths_df, bear_creek])


# crop dem and streams
for _,row in paths_df.iterrows():
    dem_file = row['dem']
    stream_ras_file = row['stream_ras']
    subbasin_file = row['subbasin']
    basinid = row['basinid']
    region = row['region']

    dem_cropped_ofile = os.path.join(odir, f"{region}_{basinid}_dem.tif")
    stream_cropped_ofile = os.path.join(odir, f"{region}_{basinid}_stream.tif")
    dem_crop = subset_raster(dem_file, subbasin_file, basinid, dem_cropped_ofile)
    stream_crop = subset_raster(stream_ras_file, subbasin_file, basinid, stream_cropped_ofile)
