import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

lipid_head_group_dict = {'popc' : 'PO4', 'pope' : 'NH3', 'dpsm' : 'NC3', 'dpg3' : 'GM1', 'chol' : 'ROH', 'pops' : 'CNO', 'pop2' : 'P1'}
lipid_group_dict = {'popc' : 'POPC', 'pope' : 'POPE', 'dpsm' : 'DPSM', 'dpg3' : 'DPG3', 'chol' : 'CHOL', 'pops' : 'POPS', 'pop2' : 'POP2'}

def main(coords, traj, lipid_sel, nbins):

    u = mda.Universe(coords, traj)

    #sel = "resname POPC and (name " + str(lipid_head_group_dict['popc'] + ")")
    sel = "resname " + str(lipid_group_dict[str(lipid_sel)]) + " and (name " + str(lipid_head_group_dict['dpg3'] + ")")

    print (sel)

    x = []
    y = []

    for frame in range(len(u.trajectory)):

        u.trajectory[frame]

        if frame % 10 == 0:

            print frame

        sel_coords = u.select_atoms(sel).positions

        L = LeafletFinder(universe=u, selectionstring=sel)

        leaflet0 = L.groups(0)
        #leaflet1 = L.groups(1)

        #for particle in range(len(sel_coords)):

        for particle in range(len(leaflet0)):

            #x.append(sel_coords[particle][0])
            #y.append(sel_coords[particle][1])

            x.append(leaflet0[particle].position[0])
            y.append(leaflet0[particle].position[1])

    xedges = np.linspace(np.min(x), np.max(x), nbins)
    yedges = np.linspace(np.min(y), np.max(y), nbins)


    H, blahx, blahy = np.histogram2d(x=x, y=y, bins=(xedges, yedges))

    plt.imshow(H, interpolation='bicubic', cmap='plasma', extent=[0, np.max(xedges), 0, np.max(yedges)])

    ax = plt.subplot(111)
    ax.set_title("GM3 Average Density (10 us)")

    ax.set_xlabel("x")
    ax.set_ylabel("y")

    plt.show()

    ##### PLOTTING #####

    # PLOT WITH 1D HISTOGRAMS

    nullfmt = NullFormatter()  # no labels

    # definitions for the axes
    # left, width = 0.1, 0.65
    # bottom, height = 0.1, 0.65
    # bottom_h = left_h = left + width + 0.02

    # rect_scatter = [left, bottom, width, height]
    # rect_histx = [left, bottom_h, width, 0.2]
    # rect_histy = [left_h, bottom, 0.2, height]
    #
    #
    #
    # # start with a rectangular Figure
    # plt.figure(1, figsize=(8, 8))
    #
    # axScatter = plt.axes(rect_scatter)
    # axHistx = plt.axes(rect_histx)
    # axHisty = plt.axes(rect_histy)
    #
    # #axHistx.xaxis.set_major_formatter(nullfmt)
    # #axHisty.yaxis.set_major_formatter(nullfmt)
    #
    # #axScatter.imshow(H, cmap='plasma', extent=[np.min(x), np.max(x), np.min(y), np.max(y)])
    #
    # axScatter.imshow(H, cmap='plasma', extent=[0, np.max(x), 0, np.max(y)])
    #
    # #axHistx.hist(x, bins=xedges)#, histtype='step')
    #
    # #axHisty.hist(y, bins=yedges, orientation='horizontal')#, histtype='step')
    #
    # axHistx.set_xlim(axScatter.get_xlim())
    # axHisty.set_ylim(axScatter.get_ylim())
    #
    #
    # #plt.colorbar()
    #
    #
    # plt.show()


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', type=str, help='System coordinate file', required=True, dest='coords')
    parser.add_argument('-x', type=str, help='System trajectory file', required=True, dest='traj')
    parser.add_argument('-l', type=str, help='Lipid selection', required=True, dest='lipid')
    args = parser.parse_args()

    main(coords=args.coords, traj=args.traj, lipid_sel=args.lipid, nbins=20)