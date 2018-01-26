import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

gro_file = "tictac_long_100ns_pbcmol_corrected_lipid_and_prot.gro"
traj_file = "pr_tictac_10msr2_skip10_pbcmol_corrected_lipid_and_prot.xtc"

u = mda.Universe(gro_file, traj_file)

sel = "name PO4"
bb_sel = "name BB"
cut_off = 5

count = 0
lipid_count = []
time_list = []

total_no_lipids = u.select_atoms(sel)

u.trajectory[0]

for frame in range(0, len(u.trajectory), 10):

    u.trajectory[frame]

    print frame


    #lipids = u.select_atoms(sel)
    #bb_atoms = u.select_atoms(bb_sel)

    lipids_pos = u.select_atoms(sel).positions
    bb_atoms_pos = u.select_atoms(bb_sel).positions

    dist = mda.analysis.distances.distance_array(bb_atoms_pos, lipids_pos, backend='OpenMP')

    p = dist[np.where(dist <= cut_off)]

    count += len(p)

    lipid_count.append(100 * (float(count) / len(total_no_lipids)))
    time_list.append(u.trajectory.time)

time_list[:] = [n / 1000000 for n in time_list]

fitx = [10, 20, 30, 40, 50]
extrap = UnivariateSpline(time_list, lipid_count, k=1)
fity = extrap(fitx)

ax = plt.subplot(111)
ax.set_xlabel("Time / microseconds", size=10)
ax.set_ylabel("Total % of Lipid Encountering Proteins", size=10)

plt.scatter(time_list, lipid_count, s=0.5, label='Current Data')
plt.plot(fitx, fity, lw=1, color='orange', linestyle='--', label='Extrapolated Data')
plt.show()
plt.savefig('test.png')