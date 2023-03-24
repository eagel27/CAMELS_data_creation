import os

SOLAR_METALLICITY = 0.0127

SIM_NAME = 'IllustrisTNG'
SIM_DIR = '1P_1_0'
BASE_PATH = '/home/jovyan/Data/Sims/{}/{}/'.format(SIM_NAME, SIM_DIR)

SNAPSHOTS = list(range(29, 34))

RESULTS_PATH = '/home/jovyan/home/'
DATASET_PATH = os.path.join(RESULTS_PATH, 'Datasets_1D')
OFFSETS_PATH = os.path.join(RESULTS_PATH, 'Offsets/{}/{}'.format(SIM_NAME, SIM_DIR))
MET_AGE_PATH = os.path.join(RESULTS_PATH, 'MetAge/{}/{}'.format(SIM_NAME, SIM_DIR))
EXSITU_PATH = os.path.join(RESULTS_PATH, 'Exsitu/{}/{}'.format(SIM_NAME, SIM_DIR))
