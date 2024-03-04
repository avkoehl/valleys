import argparse

from pyvalleys.workflow import valley_floors

parser = argparse.ArgumentParser(description="Extract valley floors")
parser.add_argument("dem", help="Path to DEM file")
parser.add_argument("flowlines", help="Path to flowlines file")
parser.add_argument("config", help="Path to config file")
parser.add_argument("wbt", help="Path to whitebox tools executable")
parser.add_argument("terrain", help="Path to terrain outputs directory")
parser.add_argument("ofile", help="Path to save valley floor shapefile")
parser.add_argument("debug", help="Path to file to store the valley wall breakpoints shapefile", default=None)
	
args = parser.parse_args()

valley_floors(args.dem, args.flowlines, args.config, args.wbt, args.terrain, args.ofile, args.debug)
