import os

import rioxarray
from shapely.geometry import Point
import whitebox

from valleys.cross_section import get_cross_section_points
from valleys.terrain import rezero_alphas, rioxarray_sample_points

def align_to_terrain(wbt, nhd_network, dem):

    # save dem to file
    dem.rio.to_raster(os.path.join(wbt.work_dir, 'dem.tif'))

    wbt.flow_accumulation_full_workflow(
            'dem.tif',
            'flow_dir.tif', 
            'flow_accum.tif',
            out_type = 'cells',
            log = False,
            clip = False,
            esri_pntr = False
            )

    # filter flow lines
    nhd_network = nhd_network.loc[nhd_network['FTYPE'] == '460']
    nhd_network = nhd_network.loc[~((nhd_network['StartFlag'] == 1) & (nhd_network['LENGTHKM'] < 1))]
    seed_points = nhd_network.loc[nhd_network['StartFlag'] == 1]['geometry'].apply(lambda x: Point(x.coords[0]))

    # save seed points to file
    seed_points.to_file(os.path.join(wbt.work_dir, 'seed_points.shp'))

    wbt.snap_pour_points(
            'seed_points.shp',
            'flow_accum.tif',
            'snap_points.shp',
            snap_dist = 50
            )

    wbt.trace_downslope_flowpaths(
            'snap_points.shp',
            'flow_dir.tif',
            'flow_paths.tif'
            )

    wbt.raster_streams_to_vector(
            'flow_paths.tif',
            'flow_dir.tif',
            'flow_paths.shp'
            )

    aligned = gpd.read_file(os.path.join(wbt.work_dir, 'flow_paths.shp'))
    return aligned


def network_cross_section_points(aligned, dem):
    # need dem to recenter the cross section points

    points_dfs = []
    for index, row in aligned.iterrows():
        flowline = row['geometry']
        centerline = flowline.simplify(100)
        points = get_cross_section_points(centerline, xs_spacing=100, xs_width=1000, xs_point_spacing=10)
        points['elevation'] = rioxarray_sample_points(dem, points)
        points = rezero_alphas(points)
        points = points.drop('elevation', axis=1)
        points['stream_id'] = row['FID']
        points_dfs.append(points)

    points = pd.concat(points_dfs, ignore_index=True)
    return points

def break_points(wbt, points, dem):
    dem.rio.to_raster(os.path.join(wbt.work_dir, 'dem.tif'))

    # curvature full workflow
    wbt.gaussian_filter("dem.tif", "dem_gauss.tif", sigma=1.5)
    wbt.slope("dem_gauss.tif", "slope.tif")
    wbt.profile_curvature("dem_gauss.tif", "profile_curvature.tif")

    points['elevation'] = rioxarray_sample_points(dem, points, 'dem.tif')
    points['slope'] = rioxarray_sample_points(dem, points, 'slope.tif')
    points['profile_curvature'] = rioxarray_sample_points(dem, points, 'profile_curvature.tif')

    break_points = get_break_points(points, peak_threshold=0.005, slope_threshold=.4)
    return break_points

def valley_polygon(wbt, break_points, aligned, dem):

    # get hand for whole area

    # get watershed for each stream
    # then for each stream, get HAND elevation threshold
    # subset hand to that threshold within the watershed
    # get polygon of that region
    # save df of stream_id, polygon
    
    # concat all the polygons together
    return valleys


# ISSUE: map nhd ids to the aligned streams
# rethink the handling of files vs rasters
