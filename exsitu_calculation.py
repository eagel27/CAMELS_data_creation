from constants import BASE_PATH

import numpy as np
import illustris_python.illustris_python as il


def get_insitu_stellar_particle_ids(basePath, subhalo_id, snapshot, previous_snapshot_time):
    """
        Calculate which particles were formed insitu for said subhalo_id at snapshot time
    """
    formed_insitu, particle_ids, masses = set(), set(), []
    fields = ['ParticleIDs', 'GFM_StellarFormationTime', 'Masses']   
    stars = il.snapshot.loadSubhalo(basePath, snapshot, subhalo_id,
                                    'stars', fields=fields)
    
    header = il.groupcat.loadHeader(basePath, snapshot)
    snapshot_time = header['Time']
    if stars['count'] > 0:
        particle_ids = stars['ParticleIDs']
        formation_time = stars['GFM_StellarFormationTime']
        masses = stars['Masses']
        formed_insitu = set(particle_ids[np.where((previous_snapshot_time < formation_time[:]) & 
                                                  (formation_time[:] <= snapshot_time))])
        particle_ids = particle_ids[formation_time[:] > 0]
        masses = masses[formation_time[:] > 0]
    return formed_insitu, snapshot_time, particle_ids, masses


def get_exsitu_fraction(subhalo_id, snapshot):
    """
        Calculate ex-situ fraction for subhalo_id at snapshot provided
    """
    fields = ['SubhaloMass','SubfindID','SnapNum']
    tree = il.sublink.loadTree(BASE_PATH, snapshot, subhalo_id, fields=fields, onlyMPB=True)

    insitu_overall = set()
    previous = 0
    for i in reversed(range(len(tree['SnapNum'][:]))):
        insitu, previous, all_ids, masses = get_insitu_stellar_particle_ids(BASE_PATH,
                                                                            tree['SubfindID'][i],
                                                                            tree['SnapNum'][i],
                                                                            previous)
        insitu_overall.update(insitu)

    insitu_overall = insitu_overall.intersection(all_ids)
    insitu_mass = sum(m for i, m in enumerate(masses) if all_ids[i] in insitu_overall)
    total_mass = sum(masses)
    
    exsitu_f = -1
    if total_mass > 0:
        exsitu_f = (total_mass - insitu_mass)/total_mass
    
    return exsitu_f

