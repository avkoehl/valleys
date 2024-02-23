import geopandas as gpd

def filter_nhd_flowlines(flowlines):
    # remove small headwater streams < 1km and non-river features
    # check for required columns
    required = ['StartFlag', 'FTYPE', 'geometry', 'LENGTHKM']

    if not all([col in nhd_network.columns for col in required]):
        raise ValueError('Input NHDPlus dataset must contain columns: StartFlag, FTYPE, geometry, LENGTHKM')

    # filter flow lines
    nhd_network = nhd_network.loc[nhd_network['FTYPE'] == 'StreamRiver']
    nhd_network = nhd_network.loc[~((nhd_network['StartFlag'] == 1) & (nhd_network['LENGTHKM'] < 1))]
    return nhd_network
