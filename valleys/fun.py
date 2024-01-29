import os

import rioxarray
from shapely.geometry import Point
import whitebox


def find_seed_points(nhd_network):
    # filter flow lines
    nhd_network = nhd_network.loc[nhd_network['FTYPE'] == '460']
    nhd_network = nhd_network.loc[~((nhd_network['StartFlag'] == 1) & (nhd_network['LENGTHKM'] < 1))]
    seed_points = nhd_network.loc[nhd_network['StartFlag'] == 1]['geometry'].apply(lambda x: Point(x.coords[0]))
    return seed_points

def align_to_terrain(wbt, nhd_network, dem_file):
    seed_points = find_seed_points(nhd_network)
    seed_points.to_file(os.path.join(wbt.work_dir, 'seed_points.shp'))

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

    paths = {
            'conditioned_dem': os.path.join(wbt.work_dir, 'conditioned_dem.tif'),
            'flow_dir': os.path.join(wbt.work_dir, 'flow_dir.tif'),
            'flow_accum': os.path.join(wbt.work_dir, 'flow_accum.tif'),
            'flow_paths': os.path.join(wbt.work_dir, 'flow_paths.tif'),
            'flow_paths_shp': os.path.join(wbt.work_dir, 'flow_paths.shp'),
            'snap_points': os.path.join(wbt.work_dir, 'snap_points.shp'),
            }
    return paths

def network_cross_section_points(aligned, xs_spacing=100, 
                                 xs_width=1000, point_spacing=10):

    points_dfs = []
    for index, row in aligned.iterrows():
        flowline = row['geometry']
        centerline = flowline.simplify(100)
        points = get_xs_points(centerline, xs_spacing, xs_width, point_spacing)
        points['stream_id'] = row['FID']
        points_dfs.append(points)

    points = pd.concat(points_dfs, ignore_index=True)
    return points

def compute_terrain_rasters(wbt, dem_file, stream_file, flow_dir_file, sigma=1.5):
    dem_file = os.path.abspath(dem_file)
    stream_file = os.path.abspath(stream_file)
    flow_dir_file = os.path.abspath(flow_dir_file)

    wbt.gaussian_filter(dem_file, "dem_gauss.tif", sigma=sigma)
    wbt.slope("dem_gauss.tif", "slope.tif")
    wbt.profile_curvature("dem_gauss.tif", "profile_curvature.tif")

    # hand
    wbt.elevation_above_stream(dem_file, stream_file, "hand.tif")

    # subbasins
    wbt.subbasins(
                flow_dir_file,
                "streams.tif" 
                "subbasins.tif",
                esri_pntr=False, 
                )

    # hillslopes
    wbt.hillslopes(
            flow_dir_file,
            "streams.tif",
            "hillslopes.tif",
            esri_pntr=False
            )

    paths = {
            'dem_gauss': os.path.join(wbt.work_dir, 'dem_gauss.tif'),
            'slope': os.path.join(wbt.work_dir, 'slope.tif'),
            'profile_curvature': os.path.join(wbt.work_dir, 'profile_curvature.tif'),
            'hand': os.path.join(wbt.work_dir, 'hand.tif'),
            'subbasins': os.path.join(wbt.work_dir, 'subbasins.tif'),
            'hillslopes': os.path.join(wbt.work_dir, 'hillslopes.tif')
            }
    return paths

def xs_boundary_points(xs_points, dem_file, curv_file):
    def boundary_point(half_profile_curvature):
        # get break points
        # filter break points to get boundary point
        return boundary_point

    # get left_boundary
    # get right_boundary
    return [{'xsid': xsid, 'left': left, 'right': right}]

#def delineate_valley(stream_id, points, dem_file, curv_file, hand_file, subbasin_file, hillslope_file):
#
#    bpoints = []
#    for xs in xsections:
#        fix_alpha(xs)
#        left, right = xs_boundary_points(xs, dem_file, curv_file)
#        bpoints.append({'xsid': xs['xsid'], 'left': left, 'right': right})
#    boundary_points = pd.concat(bpoints_dfs, ignore_index=True)
#
#    # for hillslope in hillslope with atleast 30 points:
#    # get mean hand
#    # if mean hand significantly different from mean hand of all hillslopes:
#    # keep top hillslope
#    # get .6 quantile
#
#    # subset hand to watershed
#
#    # threshold hand
#
#    # threshold slope
#
#    # vectorize + cleanup
#
#    return valley
