import os
import numpy as np
import h5py
import illustris_python.illustris_python as il

from constants import *


def get_offsets(snapshot):
    header = il.groupcat.loadHeader(BASE_PATH, snapshot)
    scaling_factor = header['Time']

    halos = il.groupcat.loadHalos(BASE_PATH, snapshot, fields=['GroupLenType', 'GroupNsubs'])
    halo_nsubs = halos['GroupNsubs']
    halo_offsets = halos['GroupLenType'].cumsum(axis=0)
    # Insert a zero-row in the beginning and remove last row
    halo_offsets = np.insert(halo_offsets, 0, 0, axis=0)[:-1, :]
    
    subhalos = il.groupcat.loadSubhalos(BASE_PATH, snapshot, fields=['SubhaloLenType'])
    subhalo_offsets = subhalos.cumsum(axis=0)
    # Insert a zero-row in the beginning and remove last row
    subhalo_offsets = np.insert(subhalo_offsets, 0, 0, axis=0)[:-1, :]
    
    # Fix offsets for the inner fuzz particles 
    halo_nsubs_cum = halo_nsubs.cumsum()
    subhalo_offsets_upd = np.copy(subhalo_offsets)
    sub_index = 0
    for hid in range(len(halo_nsubs) - 1):
        subhalos_halo = halo_nsubs[hid]
        subhalo_particles = np.sum(subhalos[sub_index:sub_index+subhalos_halo], axis=0)
        inner_fuzz = halo_offsets[hid+1] - halo_offsets[hid] - subhalo_particles
        sub_index += subhalos_halo
        
        # Update all subhalo offsets after that the last subhalo of that halo
        subhalo_offsets_upd[halo_nsubs_cum[hid]:, ] += inner_fuzz
        
    return halo_offsets, subhalo_offsets_upd


def save_offsets(snapshot, halo_offsets, subhalo_offsets):
    offset_path = os.path.join(OFFSETS_PATH, 'offsets_%03d.hdf5' % snapshot)
    with h5py.File(offset_path, 'w') as f:
        g_grp = f.create_group("Group")
        g_grp.create_dataset('SnapByType', data=halo_offsets)
        s_grp = f.create_group("Subhalo")
        s_grp.create_dataset('SnapByType', data=subhalo_offsets)
        

def create_offsets(snapshot):
    halo_offsets, subhalo_offsets = get_offsets(snapshot)
    save_offsets(snapshot, halo_offsets, subhalo_offsets)

