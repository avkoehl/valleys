import os

import rioxarray
from shapely.geometry import Point, box
from scipy.ndimage import gaussian_filter
import whitebox

from valleys.utils import chomp_raster

def clip_to_subbasin(dem, hillslopes, flow_dir, stream, subbasins, subbasin_id):
    dem_clipped = dem.where(subbasins == subbasin_id)
    dem_clipped = chomp_raster(dem_clipped)
    
    hillslope_clipped = hillslopes.where(subbasins == subbasin_id)
    hillslope_clipped = chomp_raster(hillslope_clipped)

    flow_dir_clipped = flow_dir.where(subbasins == subbasin_id)
    flow_dir_clipped = chomp_raster(flow_dir_clipped)

    stream_clipped = stream.where(stream == subbasin_id)
    stream_clipped = stream_clipped.rio.clip([box(*dem_clipped.rio.bounds())], drop=True)
    return dem_clipped, hillslope_clipped, flow_dir_clipped, stream_clipped

def flow_accumulation_workflow(wbt, dem_file):
    wbt.flow_accumulation_full_workflow(
            os.path.abspath(dem_file),
            'conditioned_dem.tif',
            'flow_dir.tif', 
            'flow_accum.tif',
            out_type = 'cells',
            log = False,
            clip = False,
            esri_pntr = False
            )
    return {
            'conditioned_dem': os.path.join(wbt.work_dir, 'conditioned_dem.tif'),
            'flow_dir': os.path.join(wbt.work_dir, 'flow_dir.tif'),
            'flow_accum': os.path.join(wbt.work_dir, 'flow_accum.tif')
            }

def find_seed_points(nhd_network):
    # filter flow lines
    nhd_network = nhd_network.loc[nhd_network['FTYPE'] == 'StreamRiver']
    nhd_network = nhd_network.loc[~((nhd_network['StartFlag'] == 1) & (nhd_network['LENGTHKM'] < 1))]
    seed_points = nhd_network.loc[nhd_network['StartFlag'] == 1]['geometry'].apply(lambda x: Point(x.coords[0]))
    return seed_points

def align_to_terrain(wbt, nhd_network, flow_acc_file, flow_dir_file):
    seed_points = find_seed_points(nhd_network)
    seed_points.to_file(os.path.join(wbt.work_dir, 'seed_points.shp'))

    wbt.snap_pour_points(
            'seed_points.shp',
            flow_acc_file,
            'snap_points.shp',
            snap_dist = 50
            )

    wbt.trace_downslope_flowpaths(
            'snap_points.shp',
            flow_dir_file,
            'flow_paths.tif'
            )

    wbt.stream_link_identifier(
            flow_dir_file,
            'flow_paths.tif',
            'flow_paths_identified.tif'
            )

    wbt.raster_streams_to_vector(
            'flow_paths_identified.tif',
            flow_dir_file,
            'flow_paths.shp'
            )

    return {
            'streams': os.path.join(wbt.work_dir, 'flow_paths_identified.tif'),
            'streams_shp': os.path.join(wbt.work_dir, 'flow_paths.shp'),
            }

def segment_subbasins(wbt, streams_file, flow_dir_file):
    # subbasins
    wbt.subbasins(
            d8_pntr = flow_dir_file,
            streams = streams_file,
            output= "subbasins.tif",
            esri_pntr=False, 
            )

    # hillslopes
    wbt.hillslopes(
            flow_dir_file,
            streams_file,
            "hillslopes.tif",
            esri_pntr=False
            )
    return {
            'subbasins': os.path.join(wbt.work_dir, 'subbasins.tif'),
            'hillslopes': os.path.join(wbt.work_dir, 'hillslopes.tif')
            }

def compute_terrain_rasters(wbt, dem_file, stream_file, flow_dir_file, sigma=1.5):
    # so that hand doesn't have gaps as occurs with breach depressions least cost
    wbt.fill_depressions(
            dem_file, # breached depressions 
            "filled_dem.tif",
            fix_flats = True,
            flat_increment = None,
            max_depth = None
            )

    # hand
    wbt.elevation_above_stream("filled_dem.tif", stream_file, "hand.tif")

    # use gaussian filter to smooth the dem, prefer scipy over wbt
    # because scipy gaussian has options to maintain the dimensions, 
    # i.e keep the same size as the input
    dem = rioxarray.open_rasterio(dem_file)
    dem.data = gaussian_filter(dem, sigma=sigma)
    dem.rio.to_raster(os.path.join(wbt.work_dir, "dem_gauss.tif"))

    #wbt.gaussian_filter(dem_file, "dem_gauss.tif", sigma=sigma)
    wbt.slope("dem_gauss.tif", "slope.tif")
    wbt.profile_curvature("dem_gauss.tif", "profile_curvature.tif")

    return {
            'dem_gauss': os.path.join(wbt.work_dir, 'dem_gauss.tif'),
            'slope': os.path.join(wbt.work_dir, 'slope.tif'),
            'curvature': os.path.join(wbt.work_dir, 'profile_curvature.tif'),
            'hand': os.path.join(wbt.work_dir, 'hand.tif'),
            }
