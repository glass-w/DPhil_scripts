import seaborn
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import distances
import matplotlib.pyplot as plt
import glob, os
from os.path import join
from datetime import datetime

#repeat = os.path.basename(os.getcwd())

#u = mda.Universe("NM_" + str(repeat) +"_resid_90_190.gro", "NM_" + str(repeat) +"_resid_90_190.xtc")

u = mda.Universe("repeat_3_analysis_resid_90_190.gro", "repeat_3_analysis_resid_90_190.xtc")

sel = "(resname ARG LYS ASP GLU) and (resid 90:190)"

resid = list(u.select_atoms(sel).atoms.residues.resids)
resnames = list(u.select_atoms(sel).resnames)

def make_hmap(universe, resid):
    ''' makes an empty numpy array to be populated
     universe = a single frame or whole trajectory
     resid = list of residue ids
     '''

    name_store = []
    resid_store = []

    for item in range(len(resid)): # Get the names and resids of the amino acids of interest

        aa = universe.select_atoms("resid " + str(resid[item]))

        if aa.resnames[0] == "ASP":

            name_store.append(aa.resnames[0])

            resid_store.append(resid[item])

        elif aa.resnames[0] == "LYS":

            name_store.append(aa.resnames[0])

            resid_store.append(resid[item])

        elif aa.resnames[0] == "GLU":

            name_store.append(aa.resnames[0])

            resid_store.append(resid[item])

        elif aa.resnames[0] == "ARG":

            name_store.append(aa.resnames[0])

            resid_store.append(resid[item])

    hmap = np.zeros((len(name_store), len(name_store))) # define the array size before making the heathmap

    return hmap, name_store, resid_store

def populate_hmap(hmap, universe, resid):
    '''populates empty heatmap array for a single frame
    hmap = empty numpy array of N x N
    universe = a single frame or whole trajectory
    resid = list of residue ids
    aa = amino acids of interest
    '''

    hmap_pop = hmap

    com_store = []

    for item in range(len(resid)):

        aa = universe.select_atoms("resid " + str(resid[item]))

        if aa.resnames[0] == "ASP":

            asp_com = aa.atoms[[8, 9]].center_of_mass()

            com_store.append(asp_com)

        elif aa.resnames[0] == "LYS":

            lys_com = aa.atoms[[16]].center_of_mass()

            com_store.append(lys_com)

        elif aa.resnames[0] == "GLU":

            glu_com = aa.atoms[[11, 12]].center_of_mass()

            com_store.append(glu_com)

        elif aa.resnames[0] == "ARG":

            arg_com = aa.atoms[[15, 16, 19]].center_of_mass()

            com_store.append(arg_com)

    contacts = distances.contact_matrix(np.array(com_store).astype(np.float32),
                                        cutoff=4.0, returntype="numpy", box=universe.dimensions)

    hmap_pop += contacts.astype(int)

    # for k in range(len(contacts)):
    #
    #     for l in range(len(contacts)):
    #
    #         if contacts[k, l] == True:
    #
    #             hmap_pop[k, l] = hmap[k, l] + 1

    hmap_pop = (hmap_pop / np.amax(hmap_pop)) * 100 # normalise to make a percentage

    return hmap_pop

def plot_hmap(data, names, resids):
    '''plot a single frame
    data = the populated hmap
    names = list of amino acid names
    resids = list of resids of amino acids
    '''

    mask = np.tri(data.shape[0], k=-1)

    data = np.ma.array(data, mask=mask)

    x = []
    y = []
    x.extend(range(len(names)))
    y.extend(range(len(names)))

    xticks = []
    yticks = []

    for i in range(len(names)):    # to get the correct x and y axis labels

        xticks.append(str(names[i]) + str(' ') + str(resids[i]))
        yticks.append(str(names[i]) + str(' ') + str(resids[i]))

    fig = plt.figure()

    ax = fig.add_subplot(111)
    #ax.set_title('Frame Number: ' + str(frame.frame))
    
    plt.xticks(range(len(names)), xticks, rotation=90)
    fig.subplots_adjust(bottom=0.2)
    plt.yticks(range(len(names)), yticks)

    plt.suptitle('Mutated Nav1.5: Residues < 4 Angstrom')
    plt.xlabel("Residue")
    plt.ylabel("Residue")

    plt.imshow(data, cmap='viridis', interpolation='spline16')

    #plt.setp(plt.plot(y, x), color='black', linewidth=6.0)
    #plt.grid(color='white',  linestyle='dotted')

    plt.grid(True)

    plt.colorbar()

    return fig

if __name__ == "__main__":

    startTime = datetime.now()

    hmap, name_store, resid_store = make_hmap(universe=u, resid=resid)

    for frame in u.trajectory:

        u.trajectory[frame.frame]

        hmap_pop = populate_hmap(hmap=hmap, universe=u, resid=resid)

        fig = plot_hmap(data=hmap_pop, names=name_store, resids=resid_store)

        #fig.savefig("frame_{:04d}.png".format((frame.frame / 10))) # need correct format for ffmpeg to read sequen'y

	fig.savefig("frame_{:04d}.png".format((frame.frame)))

    fig.clf()

    os.system('ffmpeg -i frame_%04d.png -filter:v "setpts=3.0*PTS" movie.mp4')

#    for f in glob.glob("*.png"):

 #       os.remove(f)

    print datetime.now() - startTime
