import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from matplotlib import ticker
import pylab
import seaborn as sns


gro_file = sys.argv[1]
xtc_file = sys.argv[2]

lipid_head_group_dict = {'popc' : 'PO4', 'pope' : 'NH3', 'dpsm' : 'NC3', 'dpg3' : 'GM1', 'chol' : 'ROH', 'pops' : 'CN0', 'pop2' : 'P1'}
lipid_group_dict = {'popc' : 'POPC', 'pope' : 'POPE', 'dpsm' : 'DPSM', 'dpg3' : 'DPG3', 'chol' : 'CHOL', 'pops' : 'POPS', 'pop2' : 'POP2'}
lipid_info_store = {'popc': {'name': 'POPC', 'z': [], 'colour' : (0.6, 0.6, 0.6)},
                          'pope': {'name': 'POPE', 'z': [], 'colour' : (0.30196, 0.68627, 0.290196)},
                          'dpsm': {'name': 'Sph', 'z': [], 'colour' : (0.85, 0.60, 0.90)},
                          'dpg3': {'name': 'GM3', 'z': [], 'colour' : (0.596078, 0.30588, 0.635294)},
                          'chol': {'name': 'CHOL', 'z': [], 'colour' : (1, 0.4980, 0)},
                          'pops': {'name': 'POPS', 'z': [], 'colour' : (0.6, 0.8, 1)},
                          'pop2': {'name': 'PIP2', 'z': [], 'colour' : (1.000, 1.000, 0.000)}}


def setup_syst(gro_file, xtc_file):

    universe = mda.Universe(gro_file, xtc_file)
    print(universe)

    return universe


def density(u, lipid_sel, lipid_heads, lipids, nbins):

    sel = "resname " + str(lipids[str(lipid_sel)]) + " and (name " + str(lipid_heads[str(lipid_sel)] + ")")

    print "Lipid selection = " + str(lipids[str(lipid_sel)])

    # xu, xl, yu, yl, zu, zl = [], [], [], [], [], []

    lip_coord_store = []

    for frame in range(len(u.trajectory)):

        u.trajectory[frame]

        if frame % 10 == 0:

            print frame

        L = LeafletFinder(universe=u, selectionstring=sel)

        leaflet_upper = L.groups(0)
        leaflet_lower = L.groups(1)

        for particle in range(len(leaflet_upper)):

            lip_coord_store.append((leaflet_upper[particle].position[0], leaflet_upper[particle].position[1], leaflet_upper[particle].position[2]))
            # yu.append(leaflet_upper[particle].position[1])
            # zu.append(leaflet_upper[particle].position[2])

        # for particle in range(len(leaflet_lower)):
        #
        #     xl.append(leaflet_lower[particle].position[0])
        #     yl.append(leaflet_lower[particle].position[1])

print(lip_coord_store)