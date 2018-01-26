import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder
import numpy as np
import sys, os
import matplotlib.pyplot as plt

gro_file = "tictac_long_100ns_pbcmol_corrected_lipid_and_prot.gro"
traj_file = "pr_tictac_10msr2_skip10_pbcmol_corrected_lipid_and_prot.xtc"

u = mda.Universe(gro_file, traj_file)

def main(universe, nbins):

    box = u.dimensions[0:3]

    box_x = box[0]
    box_y = box[1]

    sel = "name BB"

    # L = LeafletFinder(universe=u, selectionstring=sel)
    #
    # leaflet0 = L.groups(0)
    # leaflet1 = L.groups(1)


    #sel_coords = u.select_atoms(sel).positions

    hx = np.zeros((nbins-1))
    hy = np.zeros((nbins - 1))
    H = np.zeros((nbins-1, nbins-1))

    #for frame in (u.trajectory[-1]):



    u.trajectory[0]

    # sel_coords = u.select_atoms(sel).positions

    #xedges = np.linspace(np.min(sel_coords[:, 0]), np.max(sel_coords[:, 0]), 500)
    #yedges = np.linspace(np.min(sel_coords[:, 1]), np.max(sel_coords[:, 1]), 500)

    x = []
    y = []

#for frame in range(len(u.trajectory)-600):

    #u.trajectory[frame]

    u.trajectory[-1]

    #print frame

    sel_coords = u.select_atoms(sel).positions

    for particle in range(len(sel_coords)):

        #if sel_coords[coord][2] >= (box[2] / 2):

        x.append(sel_coords[particle][0])
        y.append(sel_coords[particle][1])

        #x.append(leaflet0[po4].position[0])
        #y.append(leaflet0[po4].position[1])

    histx = np.histogram(x, bins=nbins)

    histy = np.histogram(y, bins=nbins)

    xedges = np.linspace(np.min(x), np.max(x), nbins)
    yedges = np.linspace(np.min(y), np.max(y), nbins)

    hist, xedges, yedges = np.histogram2d(x=x, y=y, bins=(xedges, yedges), normed=False)

    # H += hist
    # hx += histx
    # hy += histy

    H = hist
    hx = histx
    hy = histy

    plt.plot(hx)


    H = H / np.max(H)

    # definitions for the axes
    # left, width = 0.1, 0.65
    # bottom, height = 0.1, 0.65
    # bottom_h = left_h = left + width + 0.02
    #
    # rect_scatter = [left, bottom, width, height]
    # rect_histx = [left, bottom_h, width, 0.2]
    # rect_histy = [left_h, bottom, 0.2, height]
    #
    # # start with a rectangular Figure
    # plt.figure(1, figsize=(8, 8))
    #
    # axScatter = plt.axes(rect_scatter)
    # axHistx = plt.axes(rect_histx)
    # axHisty = plt.axes(rect_histy)
    #
    # axScatter.imshow(H, cmap='plasma', interpolation='bicubic')







    #plt.imshow(H, cmap='plasma', interpolation='bicubic')

    # axHistx.hist(hx, bins=nbins)
    # axHisty.hist(hy, bins=nbins, orientation='horizontal')

    #ax = plt.subplot(132, title='pcolormesh: actual edges',aspect='equal')
   # X, Y = np.meshgrid(xedges, yedges)

    #ax.pcolormesh(X, Y, H)

    #plt.colorbar()

    #plt.grid(True)

    plt.show()


if __name__ == "__main__":

    main(universe=u, nbins=100)