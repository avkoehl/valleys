"""
Function to get centerline of a valley polygon 

Ensures that the start and end of the centerline are near the start and end of the stream line

1. get points equally dispersed on polygon border
2. create voronoi diagram for those points, clip to interior of polygon
3. find path along that diagram that is longest and most straight

https://observablehq.com/@veltman/centerline-labeling

Alternatively could try skeletonization or the method introduced in 
https://esurf.copernicus.org/articles/10/437/2022/
"""
import itertools

import pandas as pd
import geopandas as gpd
import networkx as nx
import numpy as np
import rioxarray
import shapely
from shapely.geometry import Point, Polygon, LineString
from shapelysmooth import taubin_smooth # prefer taubin unless need to preserve nodes
from shapelysmooth import chaikin_smooth

from pyvalleys.gis import rioxarray_sample_points

def get_flowline_centerline(polygon, flowline, flow_acc, num_points=800, smooth=True, snapped=False):
    start, end = get_inlet_and_outlet_points(flowline, flow_acc)

    # check distance between start,end and polygon boundary
    voronoi_graph = _get_interior_voronoi_network(polygon, num_points=num_points) 
    bn = get_boundary_nodes(voronoi_graph)
    
    bn['distance_to_inlet'] = bn.distance(start)
    bn['distance_to_outlet'] = bn.distance(end)
    
    sources = bn.sort_values(by='distance_to_inlet').iloc[0:10]
    targets = bn.sort_values(by='distance_to_outlet').iloc[0:10]

    if snapped:
        # ensure path starts 
        sources = sources.iloc[[0]] 
        targets = targets.iloc[[0]]

    path = get_best_path(voronoi_graph, sources, targets)
    if smooth:
        centerline = chaikin_smooth(taubin_smooth(path))
    else:
        centerline = path
    return centerline

def get_centerline(polygon, num_points, smooth=True):
    """ generic method for getting centerline of a polygon """
    voronoi_graph = _get_interior_voronoi_network(polygon, num_points=num_points)
    bn = get_boundary_nodes(voronoi_graph)
    path = get_best_path(voronoi_graph, bn, bn)

    if smooth:
        return chaikin_smooth(taubin_smooth(path))
    else:
        return path
    
def _get_interior_voronoi_network(polygon, num_points):
    # points on boundary
    points = create_points_along_boundary(polygon, num_points)
    simple = Polygon(points)

    # get voronoi 
    voronoi = shapely.voronoi_polygons(simple, only_edges=True)
    clipped = gpd.GeoSeries(voronoi).clip(simple)

    # convert to networkx graph
    lines = clipped.explode(index_parts=False).reset_index(drop=True)
    lines = gpd.GeoDataFrame(geometry=lines)
    G = lines_to_graph(lines)
    subgraphs = nx.connected_components(G)
    largest = max(subgraphs, key=len)
    G = G.subgraph(largest)
    return G

def sinuosity(linestring):
    start = linestring.interpolate(0)
    end = linestring.interpolate(1, normalized=True) 
    return linestring.length/start.distance(end)

def lines_to_graph(gdf):
    nodes = {}
    count = 0
    G = nx.Graph()
    for i,linestring in enumerate(gdf.geometry):
        c1 = linestring.coords[0]
        c2 = linestring.coords[-1]
        if c1 not in nodes:
            nodes[c1] = count
            count += 1 
        if c2 not in nodes:
            nodes[c2] = count
            count += 1
        G.add_node(nodes[c1], linestring=i, coords=c1)
        G.add_node(nodes[c2], linestring=i, coords=c2)
        G.add_edge(nodes[c1], nodes[c2], linestring=i)
    return G
    
def create_points_along_boundary(polygon, num_points): 
    boundary = polygon.boundary
    boundary_length = boundary.length
    interval_length = boundary_length / num_points
    points = []

    distance = 0
    while distance <= boundary_length:
        point = boundary.interpolate(distance)
        points.append(point)
        distance += interval_length

    points.append(points[0])
    return points

def get_boundary_nodes(g):
    boundary_nodes = [i for i in g.nodes() if len(list(g.neighbors(i))) == 1]
    records = []
    for node in boundary_nodes:
        data = g.nodes[node]
        record = {}
        record['linestring_id'] = data['linestring']
        record['node_id'] = node
        record['geometry'] = Point(data['coords'])
        records.append(record)
    bn = pd.DataFrame.from_records(records)
    bn = gpd.GeoDataFrame(bn, geometry='geometry', crs=3310)
    return bn

def get_best_path(g, sources, targets):
    combinations = list(itertools.product(sources['node_id'], targets['node_id']))
    all_paths = []
    for c in combinations:
        paths = nx.all_simple_paths(g, source=c[0], target=c[1])
        for p in paths:
            path = LineString([Point(g.nodes[n]['coords']) for n in p])
            all_paths.append(path)
    all_paths = gpd.GeoDataFrame(geometry=all_paths, crs=3310)
    
    all_paths['length'] = all_paths.length
    # sort by length
    all_paths = all_paths.sort_values(by='length', ascending=False)
    longest_k = all_paths.iloc[0:5]
    
    # get 'sinuousity' of each path (path.length / distance)
    sins = longest_k.geometry.apply(sinuosity)
    best = longest_k.loc[sins.idxmin()]
    return best.geometry

def get_inlet_and_outlet_points(flowline, flow_acc):
    coords = gpd.GeoSeries([Point(flowline.coords[0]), Point(flowline.coords[-1])]) 
    fa = rioxarray_sample_points(flow_acc, coords)
    coords = coords.iloc[fa.argsort()]
    start = coords.iloc[0]
    end = coords.iloc[1]
    return start,end
