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

import geopandas as gpd
import rioxarray
from scipy.ndimage import gaussian_filter
from shapely.geometry import box

class Watershed:
    def __init__(self, dem, flowlines, working_dir):
        self.dem = dem
        self.flowlines = flowlines

        self.working_dir = os.path.abspath(working_dir)
        self.files = {
                'dem': self.working_dir + 'dem.tif',
                'flowlines': self.working_dir + 'flowlines.shp',
                }


        self.dataset = None

        self.dem.to_raster(self.files['dem'])
        self.flowlines.to_file(self.files['flowlines'])

    def flow_accumulation_workflow(self, wbt):
        # run wbt.flow_accumulation_full_workflow
        # check to see if it worked
        # load in the outputs
        wbt.flow_accumulation_full_workflow(
                self.files['dem'],
                self.working_dir + 'conditioned_dem.tif',
                self.working_dir + 'flow_dir.tif',
                self.working_dir + 'flow_acc.tif',
                out_type = 'cells',
                log = False,
                clip = False,
                esri_pntr = False,
                )
        # check that the files are there
        if os.path.exists(self.working_dir + 'conditioned_dem.tif') and os.path.exists(self.working_dir + 'flow_dir.tif') and os.path.exists(self.working_dir + 'flow_acc.tif'):
            # load them in
            self.files['conditioned_dem'] = self.working_dir + 'conditioned_dem.tif'
            self.files['flow_dir'] = self.working_dir + 'flow_dir.tif'
            self.files['flow_acc'] = self.working_dir + 'flow_acc.tif'
            _update_dataset()

        else:
            RuntimeError('flow accumulation workflow failed')

    def align_flowlines_to_dem(self):
        files = {
                'seed_points': self.working_dir + 'seed_points.shp',
                'snapped_seed_points': self.working_dir + 'snapped_seed_points.shp',
                'flowpaths': self.working_dir + 'flowpaths.tif',
                'flowpaths_identified': self.working_dir + 'flowpaths_identified.tif',
                'flowpaths_shp': self.working_dir + 'flowpaths.shp',
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
            _update_dataset()
        else:
            RuntimeError('align flowlines to dem failed')

    def watershed_deliniation(self):
        # check if streams_file and flow_dir exist
        if not (os.path.exists(self.files['flowpaths_identified']) and os.path.exists(self.files['flow_dir'])):
            RuntimeError('Need to align flowlines to dem')

        files = {
                'subbasins': self.working_dir + 'subbasins.tif',
                'hillslopes': self.working_dir + 'hillslopes.tif'
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
            _update_dataset()
        else:
            RuntimeError('watershed delineation failed')

    def derive_terrain_attributes(self, sigma=1.2):
        files = {
                'dem_gauss': self.working_dir + 'dem_gauss.tif',
                'slope': self.working_dir + 'slope.tif',
                'curvature': self.working_dir + 'curvature.tif',
                }

        dem = rioxarray.open_rasterio(self.files['dem'])
        dem.data = gaussian_filter(dem.data, sigma = sigma)
        dem.rio.to_raster(files['dem_gauss'])

        wbt.slope(
                files['dem_gauss'],
                files['slope'],
                )
        wbt.curvature(
                files['dem_gauss'],
                files['curvature'],
                )

        exists = [os.path.exists(file) for file in files.values()]
        if all(exists):
            for key in files:
                self.files[key] = files[key]
            _update_dataset()
        else:
            RuntimeError('derive terrain attributes failed')

    def _update_dataset(self):
        self.dataset = xr.Dataset()
        for key in self.files:
            self.dataset[key] = rioxarray.open_rasterio(self.files[key], masked=True).squeeze()

    def clip_to_subbasin(self, subbasin_id):
        # check if subbasins in files
        if not os.path.exists(self.files['subbasins']):
            RuntimeError('Need to delineate watershed first')

        clipped = xr.Dataset()
        dem = self.dataset['dem']
        dem_clipped = dem.where(self.dataset['subbasins'] == subbasin_id, drop=True)
        dem_clipped = __chomp_raster(dem_clipped)
        dem_bounds = box(*dem_clipped.rio.bounds())
        clipped['dem'] = dem_clipped

        for file in self.files:
            if file != 'dem':
                raster = self.dataset[file]
                raster_clipped = raster.rio.clip(dem_bounds)
                raster_clipped = raster_clipped.rio.clip(dem_bounds, drop=True)
                clipped[file] = raster_clipped
        return clipped

    def compute_hand(self):
        # confirm dem, flowpaths, and subbasins exist
        if not (os.path.exists(self.files['dem']) and os.path.exists(self.files['flowpaths']) and os.path.exists(self.files['subbasins'])):
            RuntimeError('Need to align flowlines to dem and delineate watershed first')

        files = {
                'filled_dem': self.working_dir + 'filled_dem.tif',
                'hand': self.working_dir + 'hand.tif',
                }

        subbasin_hand_rasters = []
        for subbasin_id in np.unique(self.dataset['subbasins'].values):
            clipped = self.clip_to_subbasin(subbasin_id)

            # save dem to temp file
            clipped['dem'].rio.to_raster(self.working_dir + f'temp_dem_{subbasin_id}.tif')
            
            # fill depressions on that temp and save to temp file
            wbt.fill_depressions(
                    self.working_dir + f'temp_dem_{subbasin_id}.tif',
                    self.working_dir + f'temp_filled_dem_{subbasin_id}.tif',
                    fix_flats = True,
                    flat_increment = None,
                    max_depth = None
                    )

            # run elevation above stream and save to temp file
            wbt.elevation_above_stream(
                    self.working_dir + f'temp_filled_dem_{subbasin_id}.tif',
                    self.files['flowpaths'],
                    self.working_dir + f'temp_hand_{subbasin_id}.tif'
                    )

            # load temp file as raster and append to subbasin_hand_rasters
            hand_raster = rioxarray.open_rasterio(self.working_dir + f'temp_hand_{subbasin_id}.tif', masked=True).squeeze()

            # remove all temp files
            os.remove(self.working_dir + f'temp_dem_{subbasin_id}.tif')
            os.remove(self.working_dir + f'temp_filled_dem_{subbasin_id}.tif')
            os.remove(self.working_dir + f'temp_hand_{subbasin_id}.tif')

        # merge subbasin_hand_rasters into single raster
        hand = merge_arrays(subbasin_hand_rasters)
        hand.rio.to_raster(files['hand'])
        self.files['hand'] = files['hand']

        _update_dataset()

    def process_watershed(self):
        self.flow_accumulation_workflow()
        self.align_flowlines_to_dem()
        self.watershed_delineation()
        self.derive_terrain_attributes()
        self.compute_hand()

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
