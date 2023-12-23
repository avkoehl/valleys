"""
To be used with cross_section.py

    - Slope Threshold
        - slope_threshold                

    - Break Point Detection
        - linear regression (OLS, LAD)
        - piecewise linear regression 
        - slope change point detection
        - curvature change point detection (wavelet)
            - scipy.signal.find_peaks_cwt
"""

def slope_threshold(y, threshold):
    """
    Returns the index of the first point that exceeds the threshold
    """
    for i in range(len(y)):
        if abs(y[i]) > threshold:
            return i
    return None


