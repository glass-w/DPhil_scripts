import seaborn
import re
import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import distance_array # faster C lib
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm as CM
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import glob, os
from os.path import join
from datetime import datetime
import argparse
import collections
import networkx as nx
from itertools import count


print ' '
print "### Modules Used ####"
print ' '
print "NumPy Version: " + np.__version__
print "MDAnalysis Version: " + mda.__version__
print "Matplotlib Version: " + matplotlib.__version__
print "NetworkX Version: " + nx.__version__
print ' '
print "#####################"
print ' '

p1_network_count_dict = collections.defaultdict(dict)  # initialise the dict to store counts
p2_network_count_dict = collections.defaultdict(dict)



# p1_network = {k: [] for k in range(len(p1_residues))}
# p2_network = {(k + len(p1_residues)): [] for k in range(len(p2_residues))}

# p1_network = {k: [] for k in range(124)}
# p2_network = {(k + 124): [] for k in range(124)}

num_protein = 3
cut_off = 10

def get_universe(coord_file, traj_file):
    u = mda.Universe(coord_file, traj_file)
    return u

def setup(universe, coord_file, traj_file):

    ### this is for proteins that are the same length, may need to change this if you want to look at mixed species

    all_prot = universe.select_atoms("protein")

    all_residues = all_prot.residues

    p1_residues = universe.select_atoms("resid 1:" + str(len(all_residues) / 3)).atoms.residues
    p2_residues = universe.select_atoms("resid " + str(len(p1_residues) + 1) + ":" + str(2*(len(p1_residues)))).residues

    p1_residues_resnames = universe.select_atoms("resid 1:" + str(len(p1_residues))).atoms.residues.resnames
    p2_residues_resnames = universe.select_atoms("resid " + str(len(p1_residues) + 1) + ":" + str(2 * (len(p1_residues)))).atoms.residues.resnames

    p1_residues_resids = universe.select_atoms("resid 1:" + str(len(p1_residues))).atoms.residues.resids
    p2_residues_resids = universe.select_atoms("resid " + str(len(p1_residues) + 1) + ":" + str(2 * (len(p1_residues)))).atoms.residues.resids

    p1_resname_and_resnum = []

    for i in range(len(p1_residues)):

        p1_resname_and_resnum.append(str(p1_residues_resnames[i]) + str(p1_residues_resids[i]))

    p2_resname_and_resnum = []

    for i in range(len(p2_residues)):

        p2_resname_and_resnum.append(str(p2_residues_resnames[i]) + str(p2_residues_resids[i]))


    # This is to store what residues are near eachother
    p1_network = {k: [] for k in range(len(p1_residues))}
    p2_network = {(k + len(p1_residues)): [] for k in range(len(p2_residues))}

    # set up the centre of geometry storage arrays
    p1_cog_store = np.zeros((len(p1_residues), 3))
    p2_cog_store = np.zeros((len(p2_residues), 3))

    for residue in range(len(p1_residues)):

        p1_cog_store[residue] = u.select_atoms("resid " + str(p1_residues_resids[residue])).center_of_geometry()
        p2_cog_store[residue] = u.select_atoms("resid " + str(p2_residues_resids[residue])).center_of_geometry()

    # store the cog's in a list of lists to access in the main part of the script
    cogs = [p1_cog_store, p2_cog_store]

    return cogs, p1_network, p2_network, p1_residues, p2_residues, p1_residues_resnames, p2_residues_resnames, p1_resname_and_resnum, p2_resname_and_resnum, p1_residues_resids, p2_residues_resids

def make_network(universe, cog_store, p1_network, p2_network, p1_residues, p2_residues, p1_residues_resnames, p2_residues_resnames, p1_resname_and_resnum, p2_resname_and_resnum, p1_residues_resids, p2_residues_resids):

    # Calculate a distance matrix between protein 1 and protein 2, i.e. access the list of lists [0] = pro1, [1] = pro2

    dist = distance_array(np.array(cog_store[0]).astype(np.float32), np.array(cog_store[1]).astype(np.float32),
                       box=universe.dimensions)

    for residue in range(len(cog_store[0])):

        for i in range(len(cog_store[1])):

            if dist[residue][i] <= cut_off:

                print "d=" + str(dist[residue][i]), str(p1_residues_resnames[residue]) + str(p1_residues[residue]), str(
                    p2_residues_resnames[i]) + str(p2_residues[i])

                p1_network[residue].append(str(p1_residues_resnames[i]) + str(p2_residues_resids[i]))  #add the residues that are within cut off of protein 1 residues
                p2_network[i + len(p1_residues)].append(str(p1_residues_resnames[residue]) + str(p1_residues_resids[residue]))


    # print "P1: ", p1_network
    # print "P2: ", p2_network

    ## Get counts

    for key in p1_network:

        p1_network_count_dict[key] = collections.Counter(p1_network[key])

    for key in p2_network:

        p2_network_count_dict[key] = collections.Counter(p2_network[key])

    G = nx.Graph() # create network graph object

    for item in p1_network.values():
        if len(item) > 0:
            G.add_nodes_from(item, protein='A')
            #print item

    for item in p2_network.values():
        if len(item) > 0:
            G.add_nodes_from(item, protein='B')
            #print item

    j = 0 #counter to make sure getting correct resid

    for item in p1_network_count_dict.values():

        if len(item) > 0:
            #print item

            for i in item:
                #print "P1 Res = " + str(p1_resname_and_resnum[j]), " is near P2 Res " + str(i)

                print str(i), item
                if str(i) in p2_resname_and_resnum:
                    weight = p1_network_count_dict[j][str(i)]
                    print weight


                    G.add_edge(str(p1_resname_and_resnum[j]), str(i), weight=weight)
        j += 1


    return G

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculates the contacts between chains")
    parser.add_argument("-c", dest="gro_file", type=str, help='the coordinate file [.gro]')
    parser.add_argument("-f", dest="xtc_file", type=str,
                        help='a corrected trajectory file, pbc artifacts removed and the protein centered')
    parser.add_argument("-d", dest="cut_off_value", type=float, help='the cut-off value')
    #parser.add_argument("-p", dest="plot_type", type=str, help="define the heatplot type, i.e residue - residue or atom - atom.")
    options = parser.parse_args()

    startTime = datetime.now()

    # u, cogs, p1_network, p2_network, p1_residues, p2_residues, p1_residues_resnames, p2_residues_resnames, p1_resname_and_resnum, p2_resname_and_resnum, p1_residues_resids, p2_residues_resids = setup(coord_file=options.gro_file, traj_file=options.xtc_file)
    #
    # make_network(u, cogs, p1_network, p2_network, p1_residues, p2_residues, p1_residues_resnames, p2_residues_resnames, p1_resname_and_resnum, p2_resname_and_resnum, p1_residues_resids, p2_residues_resids)

    u = get_universe(coord_file=options.gro_file, traj_file=options.xtc_file)

    cogs, p1_network, p2_network, p1_residues, p2_residues, p1_residues_resnames, p2_residues_resnames, \
    p1_resname_and_resnum, p2_resname_and_resnum, p1_residues_resids, p2_residues_resids = setup(universe=u,
                                                                                                 coord_file=options.gro_file,
                                                                                                 traj_file=options.xtc_file)

    for frame in u.trajectory:

        u.trajectory[frame.frame]

        G = make_network(u, cogs, p1_network, p2_network, p1_residues, p2_residues, p1_residues_resnames,
                     p2_residues_resnames, p1_resname_and_resnum, p2_resname_and_resnum, p1_residues_resids,
                     p2_residues_resids)

        if u.trajectory.frame % 10 == 0:
            print "Frame = " + str(u.trajectory.frame)

        if u.trajectory.frame == 20:
            break

    # plotting

    groups = set(nx.get_node_attributes(G, 'protein').values())
    mapping = dict(zip(sorted(groups), count()))
    nodes = G.nodes()
    colors = [mapping[G.node[n]['protein']] for n in nodes]
    #pos = nx.spring_layout(G)

    print colors

    node_posA = []
    node_posB = []
    chain = nx.get_node_attributes(G, 'protein')

    ### shift positions of nodes and labels ###

    for key in chain:
        if chain[key] == 'A':
            node_posA.append(chain[key])
        elif chain[key] == 'B':
            node_posB.append(chain[key])
        else:
            continue

    i, j = 0, 0

    y1 = np.linspace(0, len(node_posA), len(node_posA))
    y2 = np.linspace(0, len(node_posB), len(node_posB))

    positions = {}
    pos_labels = {}
    x_off = 0.2

    for n, d in G.nodes(data=True):

        if d['protein'] == 'A':
            positions[n] = np.array([0, y1[i]])
            pos_labels[n] = np.array([0 - x_off, y1[i]])
            #print pos_labels
            i += 1

        elif d['protein'] == 'B':
            positions[n] = np.array([3, y2[j]])
            pos_labels[n] = np.array([3 + x_off, y2[j]])
            #print pos_labels
            j += 1

    # colors2 = range(G.number_of_edges())
    # print colors2

    #ec = nx.draw_networkx_edges(G, pos, alpha=0.2)
    #ec = nx.draw_networkx_edges(G, pos=positions, alpha=0.2, edge_color=colors2, cmap=plt.cm.Blues)
    #nc = nx.draw_networkx_nodes(G, pos=positions, alpha=0.5, nodelist=nodes, node_color=colors, node_size=300, with_labels=False)

    nc = nx.draw(G, pos=positions, nodelist=nodes, node_color=colors, node_size=300, cmap=plt.cm.jet, alpha=0.75, with_labels=False)
    nx.draw_networkx_labels(G, pos_labels)

    # plt.colorbar(nc)
    plt.axis('off')
    plt.show()

    print "Time taken = " + str(datetime.now() - startTime)