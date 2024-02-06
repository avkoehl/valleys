"""
watershed class
# code for delineating watershed and computing terrain attributes
 - dem
 - flowlines
 ----
- dataset: ['dem', 'conditioned_dem', 'smoothed_dem', 'slope', 'curvature', 'streams', 'flow_dir', 'flow_acc', 'subbasin', 'hillslopes', 'hand']

methods:
 - flow_accumulation_workflow
 - align_flowlines_to_dem 
 - watershed_delineation
 - derive_terrain_attributes
 - clip_to_subbasin (returns clipped dataset)
 - compute hand (loops through subbasins and computes hand and then merges them)
 - process_watershed (run all methods)
 """
import os
import shutil

import geopandas as gpd
import numpy as np
import rioxarray
import xarray as xr
from scipy.ndimage import gaussian_filter
from shapely.geometry import box, Point

class Watershed:
    def __init__(self, dem, flowlines, working_dir):
        self.dem = dem
        self.flowlines = flowlines

        self.working_dir = os.path.abspath(working_dir)
        if os.path.isdir(self.working_dir):
            shutil.rmtree(self.working_dir)
        os.makedirs(self.working_dir)

        self.files = {
                'dem': os.path.join(self.working_dir, 'dem.tif'),
                'flowlines': os.path.join(self.working_dir, 'flowlines.shp'),
                }


        self.dataset = None

        self.dem.rio.to_raster(self.files['dem'])
        self.flowlines.to_file(self.files['flowlines'])

    def flow_accumulation_workflow(self, wbt):
        # run wbt.flow_accumulation_full_workflow
        # check to see if it worked
        # load in the outputs
        wbt.flow_accumulation_full_workflow(
                self.files['dem'],
                os.path.join(self.working_dir, 'conditioned_dem.tif'),
                os.path.join(self.working_dir, 'flow_dir.tif'),
                os.path.join(self.working_dir, 'flow_acc.tif'),
                out_type = 'cells',
                log = False,
                clip = False,
                esri_pntr = False,
                )
        # check that the files are there
        if os.path.exists(os.path.join(self.working_dir, 'conditioned_dem.tif')) and os.path.exists(os.path.join(self.working_dir, 'flow_dir.tif')) and os.path.exists(os.path.join(self.working_dir, 'flow_acc.tif')):
            # load them in
            self.files['conditioned_dem'] = os.path.join(self.working_dir, 'conditioned_dem.tif')
            self.files['flow_dir'] = os.path.join(self.working_dir, 'flow_dir.tif')
            self.files['flow_acc'] = os.path.join(self.working_dir, 'flow_acc.tif')
            self._update_dataset()

        else:
            RuntimeError('flow accumulation workflow failed')

    def align_flowlines_to_dem(self, wbt):
        files = {
                'seed_points': os.path.join(self.working_dir, 'seed_points.shp'),
                'snapped_seed_points': os.path.join(self.working_dir, 'snapped_seed_points.shp'),
                'flowpaths': os.path.join(self.working_dir, 'flowpaths.tif'),
                'flowpaths_identified': os.path.join(self.working_dir, 'flowpaths_identified.tif'),
                'flowpaths_shp': os.path.join(self.working_dir, 'flowpaths.shp'),
                }

        seed_points = _find_seed_points(self.flowlines)
        seed_points.to_file(files['seed_points'])

        wbt.snap_pour_points(
                files['seed_points'],
                self.files['flow_acc'],
                files['snapped_seed_points'],
                snap_dist = 50)

        wbt.trace_downslope_flowpaths(
                files['snapped_seed_points'],
                self.files['flow_dir'],
                files['flowpaths'],
                )

        wbt.stream_link_identifier(
                self.files['flow_dir'],
                files['flowpaths'],
                files['flowpaths_identified'],
                )
        wbt.raster_streams_to_vector(
                files['flowpaths_identified'],
                self.files['flow_dir'],
                files['flowpaths_shp'],
                )

        # check that files exist
        exists = [os.path.exists(file) for file in files.values()]
        if all(exists):
            for key in files:
                self.files[key] = files[key]
            self._update_dataset()
            self.flowlines = gpd.read_file(files['flowpaths_shp'])
        else:
            RuntimeError('align flowlines to dem failed')

    def watershed_delineation(self, wbt):
        # check if streams_file and flow_dir exist
        if not (os.path.exists(self.files['flowpaths_identified']) and os.path.exists(self.files['flow_dir'])):
            RuntimeError('Need to align flowlines to dem')

        files = {
                'subbasins': os.path.join(self.working_dir, 'subbasins.tif'),
                'hillslopes': os.path.join(self.working_dir, 'hillslopes.tif')
                }

        wbt.subbasins(
                self.files['flow_dir'],
                self.files['flowpaths_identified'],
                files['subbasins'],
                esri_pntr = False
                )
        wbt.hillslopes(
                self.files['flow_dir'],
                self.files['flowpaths_identified'],
                files['hillslopes'],
                esri_pntr = False
                )

        exists = [os.path.exists(file) for file in files.values()]
        if all(exists):
            for key in files:
                self.files[key] = files[key]
            self._update_dataset()
        else:
            RuntimeError('watershed delineation failed')

    def derive_terrain_attributes(self, wbt, sigma=1.2):
        files = {
                'dem_gauss': os.path.join(self.working_dir, 'dem_gauss.tif'),
                'slope': os.path.join(self.working_dir, 'slope.tif'),
                'curvature': os.path.join(self.working_dir, 'curvature.tif'),
                }

        dem = rioxarray.open_rasterio(self.files['dem'])
        dem.data = gaussian_filter(dem.data, sigma = sigma)
        dem.rio.to_raster(files['dem_gauss'])

        wbt.slope(
                files['dem_gauss'],
                files['slope'],
                )
        wbt.profile_curvature(
                files['dem_gauss'],
                files['curvature'],
                )

        exists = [os.path.exists(file) for file in files.values()]
        if all(exists):
            for key in files:
                self.files[key] = files[key]
            self._update_dataset()
        else:
            RuntimeError('derive terrain attributes failed')

    def _update_dataset(self):
        self.dataset = xr.Dataset()
        for key in self.files:
            if self.files[key].endswith('.shp'):
                continue
            self.dataset[key] = rioxarray.open_rasterio(self.files[key], masked=True).squeeze()

    def clip_to_subbasin(self, subbasin_id):
        # check if subbasins in files
        if not os.path.exists(self.files['subbasins']):
            RuntimeError('Need to delineate watershed first')

        clipped = xr.Dataset()
        dem = self.dataset['dem']
        dem_clipped = dem.where(self.dataset['subbasins'] == subbasin_id, drop=True)
        dem_clipped = _chomp_raster(dem_clipped)
        dem_bounds = box(*dem_clipped.rio.bounds())
        clipped['dem'] = dem_clipped

        for data_layer in self.dataset.data_vars:
            if data_layer != 'dem':
                raster = self.dataset[data_layer]
                raster_clipped = raster.where(self.dataset['subbasins'] == subbasin_id, drop=True)  
                raster_clipped = raster_clipped.rio.clip([dem_bounds], drop=True)
                clipped[data_layer] = raster_clipped
        
        return clipped, flowlines.loc[flowlines['strm_val'] == subbasin_id]

    def compute_hand(self, wbt):
        # confirm dem, flowpaths, and subbasins exist
        if not (os.path.exists(self.files['dem']) and os.path.exists(self.files['flowpaths']) and os.path.exists(self.files['subbasins'])):
            RuntimeError('Need to align flowlines to dem and delineate watershed first')

        files = {
                'filled_dem': os.path.join(self.working_dir, 'filled_dem.tif'),
                'hand': os.path.join(self.working_dir, 'hand.tif'),
                }

        wbt.fill_depressions(
                self.files['conditioned_dem'],
                files['filled_dem'],
                fix_flats = True,
                flat_increment = None,
                max_depth = None
            )

        wbt.elevation_above_stream(
                files['filled_dem'],
                self.files['flowpaths'],
                files['hand']
                )

        if not os.path.exists(files['filled_dem']) or not os.path.exists(files['hand']):
            RuntimeError('fill depressions failed, check conditioned_dem')

        self.files['hand'] = files['hand']
        self._update_dataset()

    def process_watershed(self, wbt):
        self.flow_accumulation_workflow(wbt)
        self.align_flowlines_to_dem(wbt)
        self.watershed_delineation(wbt)
        self.derive_terrain_attributes(wbt)
        self.compute_hand(wbt)

def _find_seed_points(nhd_network):
    # filter flow lines
    nhd_network = nhd_network.loc[nhd_network['FTYPE'] == 'StreamRiver']
    nhd_network = nhd_network.loc[~((nhd_network['StartFlag'] == 1) & (nhd_network['LENGTHKM'] < 1))]
    seed_points = nhd_network.loc[nhd_network['StartFlag'] == 1]['geometry'].apply(lambda x: Point(x.coords[0]))
    return seed_points

def _chomp_raster(raster):
    raster = raster.dropna(dim='x', how='all')
    raster = raster.dropna(dim='y', how='all')
    return raster
