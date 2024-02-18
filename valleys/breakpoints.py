import geopandas as gpd
import numpy as np
import pandas as pd
from scipy import signal
from scipy.interpolate import UnivariateSpline

def filter_after_max(profile):
    # check if the profile has a peak
    # sometimes that wont be caught in filter_ridge_crossing
    # keep only points before that peak
    temp = profile.copy()
    temp['diff'] = temp['alpha'].abs()
    # sort on diff column
    temp = temp.sort_values(by='diff', ascending=True)

    # max elevation ind
    max_ind = temp['elevation'].argmax()
    keep = temp.iloc[:max_ind]
    return profile.loc[profile['point_id'].isin(keep['point_id'])]

def filter_ridge_crossing(profile):
    # check if the profile crosses a ridgeline
    # if a point has positive curvature, and there is a big jump in HAND relative to the change in elevation
    # then it is likely that the profile has crossed a ridgeline
    # remove the points after the ridgeline
    # sometimes seems that curvature is not always positive when crossing a ridgeline
    temp = profile.copy()
    temp['diff'] = temp['alpha'].abs()
    # sort on diff column
    temp = temp.sort_values(by='diff', ascending=True)

    temp['hand_diff'] = temp['hand'].diff()
    temp['elev_diff'] = temp['elevation'].diff()
    temp['ratio'] = temp['hand_diff'] / temp['elev_diff']
    temp['ratio'] = temp['ratio'].abs()

    # find the first point where the ratio is greater than 2
    for i in range(len(temp)):
        if temp['ratio'].iloc[i] > 3 and temp['hand_diff'].iloc[i] > 15:
            #            if temp['curvature'].iloc[i] > 0:
            keep = temp.iloc[:i]
            return profile.loc[profile['point_id'].isin(keep['point_id'])]
    return profile

def prepare_xsection(xs_points):
    """
    This function prepares the cross section points for analysis
    0. make sure that the cross section has the correct columns
    1. remove any duplicate points
    2. make sure that there are enough points
    3. recenter the alpha values so that the channel is at 0
    4. make sure there are points on both sides of the channel
    5. remove points that are passed a ridgeline
    6. double check that there are enough points

    The cross section points are sampled along a line perpendicular to the 
    centerline, which is generally a smoothed version of the stream centerline
    as a result, the centerpoint of the cross section may not fall on the stream

    additionally, the cross section may not have enough points to be accurately
    useful in the peak detection analysis

    the cross section may also have duplicated points if for example the point spacing is smaller than
    the dem resolution and method='nearest' was used for the sampling

    recenter alpha so that the channel is at 0
      this is tricky because the cross section may intersect with the channel at several points
      alternatively, if the spacing of points is large, then the channel may be skipped entirely
      in which case we want to recenter to the lowest point
      in either case we want to recenter to the lowest point that ALSO is closest to the middle of the 
      cross section. Consider a case where the cross section is so wide that it goes into another valley
    """
    temp = xs_points.copy()

    # 0. make sure that the cross section has the correct columns
    req = ['alpha', 'curvature', 'hand', 'point_id', 'cross_section_id', 'hillslopes', 'slope', 'strm_val']
    if not all([col in xs_points.columns for col in req]):
        raise ValueError(f'Columns missing from cross section points. Required: {req}')

    # 1. remove any duplicate points
    #  duplicated points can occur if the point spacing is smaller than the dem resolution
    #  alpha, point_id, and point will all be unique
    duplicted = temp.duplicated(subset=['elevation', 'slope', 'curvature', 'hand', 'hillslopes', 'strm_val'])
    temp = temp.loc[~duplicted]

    # 2. make sure that there are enough points
    if len(temp) < 5:
        return None

    # 3. recenter the alpha values so that the channel is at 0
    #  Aprroach
    #   Check if multiple stream values are present in the cross section
    #          If yes, something wrong with the stream clipping step

    #   Case 1: No point has the stream value
    #        Recenter the alpha values to the point with the min(hand) value
    #        This isn't ideal, but its a good start
    #   Case 2: Only one point has the stream value
    #        Recenter the alpha values to that point
    #   Case 3: Multiple points have the stream value
    #         Find the point with the min(abs(alpha)) and recenter the alpha values to that point
    
    if len(temp['strm_val'].dropna().unique()) > 1:
        raise ValueError('Cross section has multiple stream values')

    if temp['strm_val'].isna().all():
        min_point = temp['hand'].idxmin()
        min_alpha = temp['alpha'].loc[min_point]
        temp['alpha'] = temp['alpha'] - min_alpha
    elif len(temp['strm_val'].dropna().unique()) == 1:
        strm_val = temp['strm_val'].dropna().unique()[0]
        min_point = temp['hand'].loc[temp['strm_val'] == strm_val].idxmin()
        min_alpha = temp['alpha'].loc[min_point]
        temp['alpha'] = temp['alpha'] - min_alpha
    elif len(temp['strm_val'].dropna().unique()) > 1:
        stream_points = temp.loc[temp['strm_val'].notna()]
        # find stream point with min abs alpha
        min_point = stream_points['alpha'].abs().idxmin()
        min_alpha = stream_points['alpha'].loc[min_point]
        temp['alpha'] = temp['alpha'] - min_alpha

    # 4. make sure there are points on both sides of the channel
    if len(temp.loc[temp['alpha'] < 0]) == 0 or len(temp.loc[temp['alpha'] > 0]) == 0:
        return None

    # 5. remove points that are passed a ridgeline
    temp_neg = filter_ridge_crossing(temp.loc[temp['alpha'] <= 0])
    temp_neg = filter_after_max(temp_neg)
    temp_pos = filter_ridge_crossing(temp.loc[temp['alpha'] >= 0])
    temp_pos = filter_after_max(temp_pos)
    temp = pd.concat([temp_neg, temp_pos])
    # remove duplicate points
    temp = temp.drop_duplicates(subset=['point_id'])

    # 6. double check that there are enough points
    if len(temp) < 5:
        return None

    if len(temp.loc[temp['alpha'] < 0]) == 0 or len(temp.loc[temp['alpha'] > 0]) == 0:
        return None
    return temp

def find_half_profile_break_point(profile, peak_ids, slope_threshold=20):
    profile = profile.copy()
    profile = profile.reindex(profile['alpha'].abs().sort_values().index)
    profile['int_index'] = np.arange(len(profile))

    # add the first and last points
    key_point_ids = peak_ids.copy()
    key_point_ids.insert(0, profile.iloc[0]['point_id'])
    key_point_ids.append(profile.iloc[-1]['point_id'])
    key_points = profile.loc[profile['point_id'].isin(key_point_ids)]
    key_points = key_points.drop_duplicates(subset=['alpha'])

    # get mean slope between each pair of points
    # foreach pair of points, get all the points between them and get mean of slope column
    slopes = []
    for ind, (label, row) in enumerate(key_points.iterrows()):
        if ind == len(key_points) - 1:
            break
        point_id_1 = row['point_id']
        point_id_2 = key_points.iloc[ind+1]['point_id']

        pos1 = profile.loc[profile['point_id'] == point_id_1]['int_index'].values[0]
        pos2 = profile.loc[profile['point_id'] == point_id_2]['int_index'].values[0]

        points = profile.iloc[pos1:pos2]
        mean_slope = points['slope'].mean()
        slopes.append((point_id_1, mean_slope))
    slopes = pd.DataFrame(slopes, columns=['point_id', 'mean_slope'])

    # get first index where slope > slope_threshold
    # else return None
    exceeds = slopes['mean_slope'] > slope_threshold
    if exceeds.any():
        return slopes.loc[exceeds].iloc[0]['point_id']
    return 

def find_xs_break_points(df, peak_threshold=0.002, slope_threshold=20):
    """ df is the points for a single cross section 
    has columns: alpha, curvature, hand, point_id, cross_section_id, hillslopes, slope, strm_val

    # where these are point_ids
    # return (pos_break_point, neg_break_point, peak_ids)
    """

    # remove duplicates, recenter alpha, and make sure there are enough points
    df = prepare_xsection(df)

    # if the cross section is not valid, return None
    if df is None:
        return (None, None, None)

    # find the peaks in the curvature profile
    # if done after splitting, may not find the peak at the first point in the series (alpha == 0)
    # wheras if series is considered as a whole, the peak might be found, so do it before splitting
    peak_inds, _ =  signal.find_peaks(-df['curvature'], height=peak_threshold)
    peak_ids = df['point_id'].iloc[peak_inds].tolist()

    pos_bp = find_half_profile_break_point(df.loc[df['alpha'] >=0], peak_ids, slope_threshold)
    neg_bp = find_half_profile_break_point(df.loc[df['alpha'] <= 0], peak_ids, slope_threshold)

    return (pos_bp, neg_bp, peak_ids)

def find_xs_break_points_alternate(df, peak_threshold=0.002, slope_threshold=20):
    # get curvature of elevation profile not profile of curvature
    # remove duplicates, recenter alpha, and make sure there are enough points
    df = prepare_xsection(df)

    # if the cross section is not valid, return None
    if df is None:
        return (None, None, None)

    y_spl = UnivariateSpline(df['alpha'],df['elevation'],s=0,k=4)
    y_spl_2d = y_spl.derivative(n=2)

    curvatures = y_spl_2d(df['alpha'])
    peak_inds, _ = signal.find_peaks(curvatures, height=0)

    peak_ids = df['point_id'].iloc[peak_inds].tolist()

    pos_bp = find_half_profile_break_point(df.loc[df['alpha'] >=0], peak_ids, slope_threshold)
    neg_bp = find_half_profile_break_point(df.loc[df['alpha'] <= 0], peak_ids, slope_threshold)

    return (pos_bp, neg_bp, peak_ids)

# find half profile break points_alternate
# get slope from elevation not slope dem
