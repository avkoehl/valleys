import geopandas as gpd

from valleys.subbasin import Subbasin

def 

def delineate_valleys(watershed,
                     tolerance=20,
                     xs_spacing=50,
                     xs_width = 500,
                     xs_point_spacing = 10,
                     quantile = 0.75,
                     buffer = 3,
                     slope_threshold = 20):

    results = []
    redo = []
    for sid in watershed.get_subbasin_ids():
        subbasin_data, flowline = watershed.clip_to_subbasin(sid)
        subbasin_data = prep_dataset(subbasin_data)

        subbasin = Subbasin(subbasin_data, flowline, sid)
        subbasin.valley_floor_by_breakpoints_full_workflow(
                tolerance=tolerance,
                xs_spacing=xs_spacing,
                xs_width = xs_width,
                xs_point_spacing = xs_point_spacing,
                quantile = quantile,
                buffer = buffer,
                slope_threshold = slope_threshold)
        poly = subbasin.valley_floor_polygon
        threshold = subbasin.hand_threshold
        if threshold is None:
            redo.append(sid)

        if threshold is not None:
            results.append((sid, poly, threshold, quantile, buffer, slope_threshold))
    hand_thresholds = [r[2] for r in results]
    mean = sum(hand_thresholds) / len(hand_thresholds)
    for sid in redo:
        subbasin_data, flowline = watershed.clip_to_subbasin(sid)
        subbasin_data = prep_dataset(subbasin_data)
        subbasin = Subbasin(subbasin_data, flowline, sid)
        subbasin.delineate_valley_floor(buffer=buffer, slope_threshold=slope_threshold, overwrite_hand=mean)
        results.append((sid, poly, threshold, quantile, buffer, slope_threshold))

    df = gpd.GeoDataFrame(results, columns=['ID', 'floor', 'HAND', 'quantile', 'buffer', 'max_slope'], geometry='floor')
    return df

def prep_dataset(dataset):
    mapping = {
        'conditioned_dem': 'elevation',
        'flowpaths_identified': 'strm_val'
    }
    # rename bands to match mapping
    dataset = dataset.rename(mapping)
    keys = ['elevation', 'slope', 'curvature', 'strm_val', 'hillslopes', 'flow_dir', 'hand']
    dataset = dataset[keys]
    return dataset
