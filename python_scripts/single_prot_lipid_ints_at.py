import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array # faster C lib
import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd

gro_file = sys.argv[1]
xtc_file = sys.argv[2]

u = mda.Universe(gro_file, xtc_file)

resid_list = np.asarray(u.select_atoms("protein and not backbone").atoms.residues.resids)
resname_list = np.asarray(u.select_atoms("protein and not backbone").atoms.residues.resnames)


atom_list = np.asarray(u.select_atoms("protein and not backbone").positions)

cog_store = np.zeros((len(resid_list), 3))
ints = np.array([0 for i in resid_list])

for frame in u.trajectory[::100]:

    #print('Time  = {0}\r'.format(frame.time / 1000),)
    print("Time = {} ns".format(frame.time / 1000))
    p_heads = u.select_atoms("name P and resname POPC")

    N = len(resid_list)  # pre-allocate to speed up calc
    M = len(p_heads)

    for residue in range(len(resid_list)):  # For each residue, calculate the c.o.g and store

        cog_store[residue] = u.select_atoms("resid " + str(resid_list[residue])).center_of_geometry()

    d = distance_array(np.array(cog_store).astype(np.float32), np.array(p_heads.positions).astype(np.float32),
                       box=u.dimensions,
                       result=np.empty(shape=(N, M), dtype=np.float64), backend='OpenMP')

    ints += np.sum(d <= 5, axis=1)

tot_counts = np.sum(ints.flatten())

ints = (ints / tot_counts) * 100

# Focus on protein residues that interact with lipid over certain % (e.g. 10)
# Get indices of residues above certain % interaction
ints_focus_indices = [i for i, v in enumerate(ints) if v > 0]

# Get associated values
ints_focus_values = [ints[i] for i in ints_focus_indices]


# Need to correct indices for plotting, +1 due to 0 indexing. Then convert ints to strings so all plotted together
ints_focus_indices = list(map(str, ([x + 1 for x in ints_focus_indices])))

ints_focus_resids = [resname_list[int(x)-1] + " " + str(x) for x in ints_focus_indices]

print(ints_focus_resids)


# Plot bar chart
fig = plt.figure()
ax = plt.subplot(111)

plt.xlabel('% of total {} contacts '.format('POPC') + r'$\leq$' + ' 5 ' + r'$\AA$', fontsize=12)
plt.ylabel(r'$\beta$' + '3 Residue', fontsize=12)

ax.barh(ints_focus_resids, ints_focus_values)

plt.tight_layout()
plt.savefig('protein_lipid_ints.svg', format='svg', dpi=300)

