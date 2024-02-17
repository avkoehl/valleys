# config looks like:

```
wbt_path = '~/opt/WBT/'
dem_path = 'data/1801010701/dem.tif'
nhd_mr_path = 'data/1801010701/flowlines_mr.shp'
output_dir = 'data/test/'

[params]
tolerance = 20 # for the simplification of the flowline
xs_spacing = 10 # spacing of cross-sections
xs_width = 500  # total width is twice this value
xs_point_spacing = 10 # spacing of points along cross-sections
quantile = 0.75 # sets the hand threshold
buffer = 2 # add to hand threshold 
slope_threshold = 14.1 # a little more than 25% slope

peak_threshold = 0.002
bp_slope_threshold = 20
```
