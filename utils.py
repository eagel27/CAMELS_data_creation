import illustris_python.illustris_python as il
import pandas as pd


def get_snapshot_values(basePath, snapshot):

    header = il.groupcat.loadHeader(basePath, snapshot)
    scaling_factor = header['Time']
    hubble_param = header['HubbleParam']  # h.

    # total matter density in units of the critical density
    omega0 = header['Omega0']
    # density parameter corresponding to the cosmological constant
    omegaLambda = header['OmegaLambda']
    
    boxsize = header['BoxSize'] * scaling_factor / hubble_param
    return scaling_factor, hubble_param, boxsize, omega0, omegaLambda


def get_catalog_data(basePath, snapshot, scaling_factor, hubble_param):

    GroupFirstSub = il.groupcat.loadHalos(basePath, snapshot, fields=['GroupFirstSub'])
    fields = ['SubhaloMassType','SubhaloSFRinRad', 'SubhaloPos', 
              'SubhaloHalfmassRadType', 'SubhaloSpin', 'SubhaloLenType']
    subhalos = il.groupcat.loadSubhalos(basePath, snapshot, fields=fields)

    subhalos['SubhaloMassStar'] = subhalos['SubhaloMassType'][:, 4] * 1e10 / hubble_param
    subhalos['SubhaloMassGas'] = subhalos['SubhaloMassType'][:, 0] * 1e10 / hubble_param
    subhalos['SubhaloMassBH'] = subhalos['SubhaloMassType'][:, 5] * 1e10 / hubble_param
    subhalos['SubhaloMassDM'] = subhalos['SubhaloMassType'][:, 1] * 1e10 / hubble_param
    subhalos['SubhaloSpin_x'] = subhalos['SubhaloSpin'][:, 0] / hubble_param
    subhalos['SubhaloSpin_y'] = subhalos['SubhaloSpin'][:, 1] / hubble_param
    subhalos['SubhaloSpin_z'] = subhalos['SubhaloSpin'][:, 2] / hubble_param
    subhalos['SubhaloPos_x'] = subhalos['SubhaloPos'][:, 0] * scaling_factor / hubble_param
    subhalos['SubhaloPos_y'] = subhalos['SubhaloPos'][:, 1] * scaling_factor / hubble_param
    subhalos['SubhaloPos_z'] = subhalos['SubhaloPos'][:, 2] * scaling_factor / hubble_param
    subhalos['SubhaloStarsLen'] = subhalos['SubhaloLenType'][:, 4]
    subhalos['SubhaloGasLen'] = subhalos['SubhaloLenType'][:, 0]

    subhalos['SubhaloHalfmassRadStar'] = subhalos['SubhaloHalfmassRadType'][:, 4] * scaling_factor / hubble_param

    del subhalos['SubhaloSpin']
    del subhalos['SubhaloMassType']
    del subhalos['SubhaloHalfmassRadType']
    del subhalos['SubhaloPos']
    del subhalos['SubhaloLenType']
    del subhalos['count']

    tng_df = pd.DataFrame(subhalos)
    tng_df['GalaxyID'] = tng_df.index
    tng_df['Central'] = tng_df['GalaxyID'].isin(GroupFirstSub)
    tng_df = tng_df[tng_df.SubhaloMassStar > 0]
    tng_df = tng_df[tng_df.SubhaloStarsLen > 100]
    return tng_df

