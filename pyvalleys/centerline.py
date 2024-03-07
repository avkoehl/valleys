"""
Module for extracting centerline of a polygon

Options:
    1) Voronoi based / Medial Axis 
            label_centerlines 
            centerline [doesn't do filtering]
            construct_centerline is sped up centerline [doesn't do filtering]
            trimesh medial_axis https://github.com/mikedh/trimesh/blob/b9f82a52d4f804220cb8829872dea9ba04fd6f82/trimesh/path/polygons.py#L371

            https://observablehq.com/@veltman/centerline-labeling

        Basically:
            Sample points along the border
            create voronoi polygons
            clip to polygon
            know that medial axis is path along this network
            find all the paths between all the nodes along the polygon border
            keep longest straightest path

    2) Topological Skeleton (skeletonize scikit-image)

    3) Fiona J Clubb et al. (2022) https://esurf.copernicus.org/articles/10/437/2022/

    4) Alternate Medial Axis
        = select path along voronoi diagram that starts nearest to
          channel start ands at point nearest to channel end
          
            

how does arc do it? https://pro.arcgis.com/en/pro-app/latest/tool-reference/topographic-production/polygon-to-centerline.htm

https://pygeoops.readthedocs.io/en/latest/



"""
from label_centerlines import get_centerline

import pygeoops

import matplotlib.pyplot as plt
from centerline.geometry import Centerline

import geopandas as gpd
import shapely

def construct_centerlines(input_geometry, interpolation_distance=0.5):
    borders = input_geometry.segmentize(interpolation_distance)
    voronoied = shapely.voronoi_polygons(borders, only_edges=True)

    centerlines = gpd.sjoin(gpd.GeoDataFrame(geometry=gpd.GeoSeries(voronoied.geoms)), gpd.GeoDataFrame(geometry=gpd.GeoSeries(input_geometry)), op="within")
    return centerlines.unary_union

def find_paths(centerline, start_points, end_points):
    # convert to network x
    # get list of paths
    # get lengths
    # keep top 5 paths within 20% of the distance of the longest path
    return paths

def select_straightest(paths):
    # run on just 5 paths
    # compute sum of all internal angles
    # keep straightest
    return path

valleys = gpd.read_file("valley_floors.shp")

valley = valleys.loc[30].to_frame().T
valley = gpd.GeoDataFrame(valley, geometry='geometry', crs=3310)

centerline = Centerline(valley.geometry.iloc[0], interpolation_distance=50)
geoms = centerline.geometry.geoms
test = gpd.GeoSeries(geoms)

fig, ax = plt.subplots()
valley.plot(ax=ax, facecolor='none', edgecolor='red')
test.plot(ax=ax)


centerline2 = get_centerline(valley.geometry.iloc[0], max_points=20)
test2 = gpd.GeoSeries(centerline2, crs=3310)

fig, ax = plt.subplots()
valley.plot(ax=ax, facecolor='none', edgecolor='red')
test2.plot(ax=ax)

centerline3 = pygeoops.centerline(valley.geometry.iloc[0])
test3 = gpd.GeoSeries(centerline3, crs=3310)

fig, ax = plt.subplots()
valley.plot(ax=ax, facecolor='none', edgecolor='red')
test3.plot(ax=ax)


centerline4 = construct_centerline(valley.geometry.iloc[0], interpolation_distance = 10000)
test4 = gpd.GeoSeries(centerline4, crs=3310)

fig, ax = plt.subplots()
valley.plot(ax=ax, facecolor='none', edgecolor='red')
test4.plot(ax=ax)
