from utils import *
from itertools import repeat
from exsitu_calculation import *
from offsets import *
from calculate_mw_met_age import  *
from constants import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from multiprocessing import Pool
import tqdm


def calc_offsets_sim():
    if not os.path.exists(OFFSETS_PATH):
        os.makedirs(OFFSETS_PATH)

    pool = Pool(processes=4)
    results = list(tqdm.tqdm(pool.starmap(create_offsets, zip(range(34))),
                             total=34))
    pool.close()
    pool.join()

    
def calc_met_age_sim():
    if not os.path.exists(MET_AGE_PATH):
        os.makedirs(MET_AGE_PATH)

    for snapshot in SNAPSHOTS:
        (scaling_factor, hubble_param,
            boxsize, omega0, omegaLambda) = get_snapshot_values(BASE_PATH, snapshot)
        tng = get_catalog_data(BASE_PATH, snapshot, scaling_factor, hubble_param)
        galaxy_ids = list(tng.GalaxyID)

        pool = Pool(processes=4)
        results = list(tqdm.tqdm(pool.starmap(get_mass_weighted_age_met, 
                                              zip(galaxy_ids, repeat(snapshot))),
                                 total=len(galaxy_ids)))

        pool.close()
        pool.join()

        met, age = map(list, zip(*results))
        met_age_dict = {'GalaxyID': galaxy_ids, 'Metallicity': met, 'Age': age}
        df = pd.DataFrame(met_age_dict)

        save_path = os.path.join(MET_AGE_PATH, 'met_age_%03d.hdf5' % snapshot)
        df.to_hdf(save_path, '/data')

        
def calc_exsitu_sim():
    if not os.path.exists(EXSITU_PATH):
        os.makedirs(EXSITU_PATH)

    for snapshot in SNAPSHOTS:
        (scaling_factor, hubble_param,
            boxsize, omega0, omegaLambda) = get_snapshot_values(BASE_PATH, snapshot)
        tng = get_catalog_data(BASE_PATH, snapshot, scaling_factor, hubble_param)
        galaxy_ids = list(tng.GalaxyID)

        pool = Pool(processes=4)
        results = list(tqdm.tqdm(pool.starmap(get_exsitu_fraction, 
                                              zip(galaxy_ids, repeat(snapshot))),
                                 total=len(galaxy_ids)))

        pool.close()
        pool.join()

        exsitu_dict = {'GalaxyID': galaxy_ids, 'ExsituF': results}
        df = pd.DataFrame(exsitu_dict)

        save_path = os.path.join(EXSITU_PATH, 'exsitu_%03d.hdf5' % snapshot)
        df.to_hdf(save_path, '/data')


if __name__ == '__main__':
    if not os.path.exists(OFFSETS_PATH):
        print('Calculating offsets for sim {}'.format(SIM_DIR))
        calc_offsets_sim()

    if not os.path.exists(EXSITU_PATH):
        print('Calculating exsitu for sim {}'.format(SIM_DIR))
        calc_exsitu_sim()

    if not os.path.exists(MET_AGE_PATH):
        print('Calculating MW Met and Age for sim {}'.format(SIM_DIR))
        calc_met_age_sim()

    all_df = pd.DataFrame()
    for snapshot in SNAPSHOTS:
        (scaling_factor, hubble_param, boxsize,
            omega0, omegaLambda) = get_snapshot_values(BASE_PATH, snapshot)
        tng_df = get_catalog_data(BASE_PATH, snapshot, scaling_factor, hubble_param)

        exsitu_path = os.path.join(EXSITU_PATH, 'exsitu_%03d.hdf5' % snapshot)
        ex_df = pd.read_hdf(exsitu_path, '/data')

        met_age_path = os.path.join(MET_AGE_PATH, 'met_age_%03d.hdf5' % snapshot)
        met_age_df = pd.read_hdf(met_age_path, '/data')

        tng_df = tng_df.merge(ex_df)
        tng_df = tng_df.merge(met_age_df)
        tng_df['Snapshot'] = snapshot
        all_df = all_df.append(tng_df, ignore_index=True)

    final_df = pd.DataFrame()
    for column in ('GalaxyID', 'Snapshot', 'ExsituF', 'Central', 'Age', 'Metallicity'):
        final_df[column] = all_df[column]

    total_spin_x = all_df.SubhaloSpin_x
    total_spin_y = all_df.SubhaloSpin_y
    total_spin_z = all_df.SubhaloSpin_z
    total_spin = np.sqrt(total_spin_x**2 + total_spin_y**2 + total_spin_z**2)

    final_df['Stellar_Mass'] = all_df['SubhaloMassStar']
    final_df['Gas_Mass'] = all_df['SubhaloMassGas']
    final_df['DM_Mass'] = all_df['SubhaloMassDM']
    final_df['BH_Mass'] = all_df['SubhaloMassBH']
    final_df['HMSR'] = all_df['SubhaloHalfmassRadStar']
    final_df['TotalSpin'] = total_spin

    plt.hist(final_df.ExsituF)
    plt.show()

    if not os.path.exists(DATASET_PATH):
        os.makedirs(DATASET_PATH)

    save_path = os.path.join(DATASET_PATH, 'dataset_{}.hdf5'.format(SIM_DIR))
    final_df.to_hdf(save_path, '/data')





