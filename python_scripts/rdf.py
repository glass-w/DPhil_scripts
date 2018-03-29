import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import spline
import numpy as np
from glob import glob
from datetime import datetime

print ' '
print "### Modules Used ####"
print ' '
print "NumPy Version: " + np.__version__
print "MDAnalysis Version: " + mda.__version__
print "Matplotlib Version: " + matplotlib.__version__
print ' '
print "#####################"
print ' '

startTime = datetime.now()

u = mda.Universe("../../input/alpha_1micros_randomise_r0.gro",
                 "../traj_and_gro_correct/nav_3x3_10us_r0_mixed_lipid_skip10.xtc")

lipid_head_group_dict = {'popc' : 'PO4', 'pope' : 'NH3', 'dpsm' : 'NC3', 'dpg3' : 'GM1', 'chol' : 'ROH', 'pops' : 'CN0', 'pop2' : 'P1'}
lipid_group_dict = {'popc' : 'POPC', 'pope' : 'POPE', 'dpsm' : 'DPSM', 'dpg3' : 'DPG3', 'chol' : 'CHOL', 'pops' : 'POPS', 'pop2' : 'POP2'}

bb = u.select_atoms("name BB or name SC*")

prot = u.select_atoms("protein")

print prot

popc = u.select_atoms("resname POPC")
pope = u.select_atoms("resname POPE")
dpsm = u.select_atoms("resname DPSM")
dpg3 = u.select_atoms("resname DPG3")
chol = u.select_atoms("resname CHOL")
pops = u.select_atoms("resname POPS")
pop2 = u.select_atoms("resname POP2")

lipids = [popc, pope, dpsm, dpg3, chol, pops, pop2]

i = 0

for lipid in lipids:

    print ("lipid " + str(i) + " of " + str(len(lipids)))

    lip_rdf = InterRDF(bb, lipid, 100, (0.0, 80.0))#, start=0, stop=2)
    lip_rdf.run()

    lip_rdf.bins = lip_rdf.bins / 10

    xnew = np.linspace(lip_rdf.bins.min(), lip_rdf.bins.max(), 500)  # No. represents number of points to make between T.min and T.max

    power_smooth = spline(lip_rdf.bins, lip_rdf.rdf, xnew)

    plt.plot(xnew, power_smooth, label=lipid.atoms.residues.resnames[0], alpha=0.75)

    i += 1

# lip_rdf = InterRDF(prot, dpsm, 100, (0.0, 80.0), start=0, stop=2)
# lip_rdf.run()
#
# lip_rdf.bins = lip_rdf.bins / 10
#
# xnew = np.linspace(lip_rdf.bins.min(), lip_rdf.bins.max(), 500)  # No. represents number of points to make between T.min and T.max
#
# power_smooth = spline(lip_rdf.bins, lip_rdf.rdf, xnew)
#
# plt.plot(xnew, power_smooth, label='dpg3', alpha=0.75)

def plot_multiple():

    files = glob("*.csv")

    ax = plt.subplot(111)
    ax.set_title(r"RDF of $\alpha$ subunit (mixed lipid)")

    # Temporary fix!
    #mean = np.mean(np.array([split_all[0], split_all[1]]), axis=0)

    ### Plotting ###

    y_list = []
    avg = []

    # plot repeats

    print files
    from numpy import genfromtxt

    for repeat in range(len(files)):

        plot_data = genfromtxt(files[repeat], delimiter=',')

        print plot_data

        time = plot_data[:, 0]
        rmsd = plot_data[:, 1]

        #ax.plot((time / 1000), rmsd, lw=1, alpha=0.25, label='Run ' + str(repeat))


        # if repeat == 0:
        #     label='TMD'
        #
        # elif repeat == 1:
        #     label='All'
        #
        # elif repeat == 2:
        #     label='ECD'

        ax.plot((time / 1000), rmsd, lw=1, alpha=0.75, label=label)

        y_list.append(rmsd)

    # plot average line

    # if plot_average == "y":
    #
    #     for row in range(len(rmsd)):
    #
    #         # take the same particle from each repeat
    #         particle_list = [item[row] for item in y_list]
    #
    #         # take the average over each run
    #         avg.append(np.average(particle_list))

        #rolling_mean = running_mean(avg, 100)

        #ax.plot(time / 1000, avg[:], lw=0.5, alpha=1, color='black', label='Average')

        #ax.plot(time / 1000, avg, lw=0.5, alpha=1, color='black', label='Average')

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("g(r)")
    ax.legend()
    plt.draw()

    ax.figure.savefig("multiple_RMSD.svg", format='svg', dpi=300)

    print "Plot complete."





plt.axhline(y=1, xmin=0.0, xmax=80.0, linestyle='dashed', color='black', alpha=0.25)
plt.legend()
plt.title("RDF of lipids around alpha subunit")
plt.xlabel("r (nm)")
plt.ylabel("g(r)")
plt.savefig('RDf_of_all_r_8nm_all_frames.svg', dpi=300, format='svg')


print "Time taken = " + str(datetime.now() - startTime)

plt.show()
