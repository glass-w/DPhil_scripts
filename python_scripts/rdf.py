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

# u = mda.Universe("../../input/alpha_1micros_randomise_r0.gro",
#                  "../traj_and_gro_correct/nav_3x3_10us_r0_mixed_lipid_skip10.xtc")
#
# lipid_head_group_dict = {'popc' : 'PO4', 'pope' : 'NH3', 'dpsm' : 'NC3', 'dpg3' : 'GM1', 'chol' : 'ROH', 'pops' : 'CN0', 'pop2' : 'P1'}
#lipid_group_dict = {'popc' : 'POPC', 'pope' : 'POPE', 'd# psm' : 'DPSM', 'dpg3' : 'DPG3', 'chol' : 'CHOL', 'pops' : 'POPS', 'pop2' : 'POP2'}
#
# bb = u.select_atoms("name BB or name SC*")
#
# prot = u.select_atoms("protein")
#
# print prot
#
# popc = u.select_atoms("resname POPC")
# pope = u.select_atoms("resname POPE")
# dpsm = u.select_atoms("resname DPSM")
# dpg3 = u.select_atoms("resname DPG3")
# chol = u.select_atoms("resname CHOL")
# pops = u.select_atoms("resname POPS")
# pop2 = u.select_atoms("resname POP2")
#
# lipids = [popc, pope, dpsm, dpg3, chol, pops, pop2]
#
# i = 0
#
# for lipid in lipids:
#
#     print ("lipid " + str(i) + " of " + str(len(lipids)))
#
#     lip_rdf = InterRDF(bb, lipid, 100, (0.0, 80.0))#, start=0, stop=2)
#     lip_rdf.run()
#
#     lip_rdf.bins = lip_rdf.bins / 10
#
#     xnew = np.linspace(lip_rdf.bins.min(), lip_rdf.bins.max(), 500)  # No. represents number of points to make between T.min and T.max
#
#     power_smooth = spline(lip_rdf.bins, lip_rdf.rdf, xnew)
#
#     plt.plot(xnew, power_smooth, label=lipid.atoms.residues.resnames[0], alpha=0.75)
#
#     i += 1





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

    data = ['dpg3', 'pip2', 'chol', 'pops', 'pope', 'popc', 'dpsm']

    lipid_group_dict = {'popc': 'POPC', 'pope': 'POPE', 'dpsm': 'DPSM', 'dpg3': 'GM3', 'chol': 'CHOL',
                        'pops': 'POPS', 'pip2': 'PIP$_2$'}

    #cb_friendly_cp = ['#40004b','#9970ab','#c2a5cf','#a6dba0','#5aae61','#1b7837','#00441b']

    cb_friendly_cp = ['#bfd3e6', '#9ebcda', '#8c96c6', '#8c6bb1', '#88419d', '#810f7c', '#4d004b']
    cpal = cb_friendly_cp[::-1]
    #print cpal

    colour = 0
    for df in data:

        #print df

        x, y = [], []

        with open(str(df) + ".xvg") as f:
            i = 0
            for line in f:
                i += 1

                if i >= 17:

                    cols = line.split()

                    #print cols

                    if len(cols) == 2:
                        if float(cols[0]) <= 8:
                            x.append(float(cols[0]))
                            y.append(float(cols[1]))

        ax = plt.subplot(111)

        ax.plot(x, y, lw=1.5, alpha=.75, label=lipid_group_dict[df], c=cpal[colour])
        plt.axhline(y=1, xmin=0.0, xmax=8.0, linestyle='dashed', color='black', alpha=0.25, lw=0.5)

        plt.draw()

        colour += 1

    ax.set_title(r"RDF of $\alpha$ - subunit", size=16)
    ax.set_xlabel("r (nm)", size=16)
    ax.set_ylabel("g(r)", size=16)
    ax.legend()
    ax.figure.savefig("RDF_Plot.svg", format='svg', dpi=300)

    print "Plot complete."

if __name__ == "__main__":

    plot_multiple()

# plt.axhline(y=1, xmin=0.0, xmax=80.0, linestyle='dashed', color='black', alpha=0.25)
# plt.legend()
# plt.title("RDF of lipids around alpha subunit")
# plt.xlabel("r (nm)")
# plt.ylabel("g(r)")
# plt.savefig('RDf_of_all_r_8nm_all_frames.svg', dpi=300, format='svg')
#
#
# print "Time taken = " + str(datetime.now() - startTime)
#
# plt.show()
