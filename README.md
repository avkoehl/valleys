workflow:

1) get stream network data
    - split into 1km segments
    - get additional data to be used in stratified sampling (e.g confinment, stream order,...)
2) create stratified sample of stream network
3) compute GCS metrics
    - need to access DEM
    - clip to 1km segment and different breaks from centerline (25 deg slope)
    - compute metrics on clipped DEM
4) cluster
