import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import scipy.stats as stats

def init(coords, traj):

    universe = mda.Universe(coords, traj)

    lipid_head_group_dict = {'popc' : 'PO4', 'pope' : 'NH3', 'dpsm' : 'NC3', 'dpg3' : 'GM1', 'chol' : 'ROH', 'pops' : 'CN0', 'pop2' : 'P1'}
    lipid_group_dict = {'popc' : 'POPC', 'pope' : 'POPE', 'dpsm' : 'DPSM', 'dpg3' : 'DPG3', 'chol' : 'CHOL', 'pops' : 'POPS', 'pop2' : 'POP2'}
    lipid_head_z_store = {'popc': {'name': 'POPC', 'z': []},
                          'pope': {'name': 'POPE', 'z': []},
                          'dpsm': {'name': 'DPSM', 'z': []},
                          'dpg3': {'name': 'DPG3', 'z': []},
                          'chol': {'name': 'CHOL', 'z': []},
                          'pops': {'name': 'POPS', 'z': []},
                          'pop2': {'name': 'POP2', 'z': []}}

    return universe, lipid_head_group_dict, lipid_group_dict, lipid_head_z_store

def density(u, lipid_sel, lipid_heads, lipids, nbins):

    sel = "resname " + str(lipids[str(lipid_sel)]) + " and (name " + str(lipid_heads[str(lipid_sel)] + ")")

    print "Lipid selection = " + str(lipids[str(lipid_sel)])

    xu, yu, xl, yl = [], [], [], []

    for frame in range(len(u.trajectory)):

        u.trajectory[frame]

        # if u.trajectory.frame == 100:
        #     break

        if frame % 10 == 0:

            print frame

        L = LeafletFinder(universe=u, selectionstring=sel)

        leaflet_upper = L.groups(0)
        leaflet_lower = L.groups(1)

        for particle in range(len(leaflet_upper)):

            xu.append(leaflet_upper[particle].position[0])
            yu.append(leaflet_upper[particle].position[1])

        for particle in range(len(leaflet_lower)):

            xl.append(leaflet_lower[particle].position[0])
            yl.append(leaflet_lower[particle].position[1])

    xuedges = np.linspace(np.min(xu), np.max(xu), nbins)
    yuedges = np.linspace(np.min(yu), np.max(yu), nbins)

    xledges = np.linspace(np.min(xl), np.max(xl), nbins)
    yledges = np.linspace(np.min(yl), np.max(yl), nbins)


    H_upper, plot_xu, plot_yu = np.histogram2d(x=xu, y=yu, bins=(xuedges, yuedges))
    H_lower, plot_xl, plot_yl = np.histogram2d(x=xl, y=yl, bins=(xledges, yledges))

    # compress data for return
    H = [H_upper, H_lower]
    edges = {'ul_xedge' : xuedges, 'ul_yedge' : yuedges, 'll_xedge' : xledges, 'll_yedge' : yledges}

    return lipid_sel, H, edges

def mbrane_density(u, lipid_heads, lipid_heads_z, lipids):

    #u.trajectory[-1]

    half_box = u.dimensions[2] / 2

    for frame in range(len(u.trajectory)):

        u.trajectory[frame]

        for lipid in lipids:
            print lipid

            sel = "resname " + str(lipids[lipid]) + " and (name " + str(lipid_heads[lipid] + ")")

            print sel

            L = LeafletFinder(universe=u, selectionstring=sel, pbc=False)

            print L

            leaflet_upper = L.groups(0)
            leaflet_lower = L.groups(1)

            com_upper = leaflet_upper.center_of_geometry()

            for particle in range(len(leaflet_upper)):
                lipid_heads_z[lipid]['z'].append(leaflet_upper[particle].position[2])

            for particle in range(len(leaflet_lower)):
                lipid_heads_z[lipid]['z'].append(leaflet_lower[particle].position[2])

    lipid_z_populated = lipid_heads_z

    return lipid_z_populated

def plot_2D(lipid_selection, lipid_dict, data, plot_edges):

    # PLOT 2D HEATMAP

    ## Plot upper leaflet

    ax1 = plt.subplot(121)

    ax1.set_title(str(lipid_dict[str(lipid_selection)]) + " Average Density (10 us) - Upper Leaflet")

    ax1.set_xlabel("x (nm)")
    ax1.set_ylabel("y (nm)")

    ax1.imshow(data[0], interpolation='bicubic', cmap='plasma', extent=[0, np.max(plot_edges['ul_xedge']), 0, np.max(plot_edges['ul_yedge'])])

    ## Plot lower leaflet

    ax2 = plt.subplot(122)

    ax2.set_title(str(lipid_dict[str(lipid_selection)]) + " Average Density (10 us) - Lower Leaflet")

    ax2.set_xlabel("x (nm)")
    ax2.set_ylabel("y (nm)")

    ax2.imshow(data[1], interpolation='bicubic', cmap='viridis', extent=[0, np.max(plot_edges['ll_xedge']), 0, np.max(plot_edges['ll_yedge'])])

    plt.savefig('2d_density_' + str(lipid_sel), dpi=300)

    plt.clf()

def plot_1D(u, z_data, nbins):

    # PLOT 1D CROSS SEC

    ax = plt.subplot(111)

    ax.set_title(" Average Density (at 10 us)")

    ax.set_xlabel("Density (counts)")
    ax.set_ylabel("z (nm)")

    if args.p1d_selection == 'all':

        for lipid in z_data: # z_data is a double dict'y

            data = z_data[lipid]['z']

            data[:] = [x / 10 for x in data]

            density = stats.gaussian_kde(data)

            x = np.linspace(0, (u.dimensions[2] / 10), nbins)

            ax.plot(density(x), x, label=lipid, alpha=0.75)
            # ax.axhline(u.dimensions[2] / 20)

            #ax.hist(z_data[lipid]['z'], histtype='stepfilled', orientation='horizontal', bins=nbins, label=lipid)

    else:

        data = z_data[args.p1d_selection]['z']

        data[:] = [x / 10 for x in data]

        density = stats.gaussian_kde(data)

        x = np.linspace(0, (u.dimensions[2] / 10), nbins)

        ax.plot(density(x), x, label=args.p1d_selection, alpha=0.75)
        # ax.axhline(u.dimensions[2] / 20)

        # ax.hist(z_data[lipid]['z'], histtype='stepfilled', orientation='horizontal', bins=nbins, label=lipid)



    ax.legend()

    plt.savefig("Densityavg_" + str(args.p1d_selection), dpi=300)

    plt.show()
    plt.clf()

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', type=str, help='System coordinate file', required=True, dest='coords')
    parser.add_argument('-x', type=str, help='System trajectory file', required=True, dest='traj')
    parser.add_argument('-l', type=str, help='Lipid selection - for example popc, pope, dpsm, chol etc',
                        required=False, dest='lipid')
    parser.add_argument('-p', type=str, help='Plot type: 2D histogram of density over all frames = "2d". 1D cut of z-axis of last frame = 1d', required=True, dest='plot')
    parser.add_argument('-p1d_sel', type=str, help='If plotting 1D plot, choose the lipid you want to plot e.g. "popc", default = "all"', default='all', dest='p1d_selection')
    args = parser.parse_args()

    universe, lipid_heads, lipids, lipid_heads_z = init(coords=args.coords, traj=args.traj)

    if args.plot == '2d':

        lipid_sel, twoD_density_data, twoD_plot_edges = density(u=universe, lipid_heads=lipid_heads, lipids=lipids, lipid_sel=args.lipid, nbins=20)

        plot_2D(lipid_selection=lipid_sel, lipid_dict=lipids, data = twoD_density_data, plot_edges=twoD_plot_edges)

    elif args.plot == '1d':

        lipid_z_data = mbrane_density(u=universe, lipid_heads=lipid_heads, lipid_heads_z=lipid_heads_z, lipids=lipids)

        plot_1D(u=universe, z_data=lipid_z_data, nbins=5000)

