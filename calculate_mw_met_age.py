import numpy as np
import illustris_python.illustris_python as il

from utils import get_snapshot_values
from constants import SOLAR_METALLICITY


def calculate_age(expansionFactor, hubble_param, omega0, omegaLambda):
    hubble_time = 9.778 / hubble_param  # Gyr

    factor1 = 2.0 / 3.0 / np.sqrt(omegaLambda)
    term1 = np.sqrt(omegaLambda / omega0) * expansionFactor ** 1.5
    term2 = np.sqrt(1 + omegaLambda / omega0 * expansionFactor ** 3)
    factor2 = np.log(term1 + term2)

    time = factor1 * factor2 * hubble_time

    return time


def get_mass_weighted_age_met(base_path, subhalo_id, snapshot):
    (scaling_factor, hubble_param,
        boxsize, omega0, omegaLambda) = get_snapshot_values(base_path, snapshot)

    mw_age, mw_met = 0, 0
    fields = ['ParticleIDs', 'GFM_StellarFormationTime', 'GFM_Metallicity', 'Masses']   
    stars = il.snapshot.loadSubhalo(base_path, snapshot, subhalo_id,
                                    'stars', fields=fields)
    
    if stars['count'] > 0:
        particle_ids = stars['ParticleIDs']
        formation_time = stars['GFM_StellarFormationTime']
        
        stars_mask = np.where(formation_time > 0)[0]
        masses = stars['Masses'][stars_mask]
        metallicity = stars['GFM_Metallicity'][stars_mask] / SOLAR_METALLICITY
        formation_time = formation_time[stars_mask]
        
        weighted_metallicity = np.sum(np.multiply(masses, metallicity))
        mw_met = np.divide(weighted_metallicity, np.sum(masses))
        
        stellar_ages = (calculate_age(scaling_factor, hubble_param, omega0, omegaLambda) - 
                        calculate_age(formation_time, hubble_param, omega0, omegaLambda))
        weighted_age = np.sum(np.multiply(masses, stellar_ages))
        mw_age = np.divide(weighted_age, np.sum(masses))

    return mw_met, mw_age






