"""
input: 
    valley_floors
    flow_dir
    flow_acc
    flowpaths

output:
    flowpaths_reaches 

steps:
    1. get centerlines
    2. get continous width measurements
    3. add points where significant width change
    4. remove points if segment too small
    5. clip points to flowpaths
    6. reidentify flowpaths
"""
