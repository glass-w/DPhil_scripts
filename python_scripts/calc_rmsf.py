import MDAnalysis as mda
from MDAnalysis.analysis import rms
import numpy as np
import matplotlib.pyplot as plt
import seaborn.apionly as sns
import argparse


def load_syst(coord, traj):

    universe = mda.Universe(coord, traj)

    return universe

def run_rmsf(universe):

    prot = universe.select_atoms("protein")
    R = rms.RMSF(prot)
    R.run()

    CA = prot.select_atoms("name CA")
    CA_rmsf = R.rmsf[np.in1d(prot.ids, CA.ids)]


    print(CA.resids)
    print(CA_rmsf)

    return CA, CA_rmsf

def plot_rmsf(data, resid_data):

    col = ['#8c96c6', '#8856a7', '#810f7c']


    plt.style.use('ggplot')
    sns.set_style('ticks')

    ax = plt.subplot(111)
    ax.plot(resid_data.resids, data, '-', linewidth=1, color=col[-1])
    ax.fill_between(resid_data.resids, data, alpha=0.1, color=col[-1])

    sns.despine(ax=ax, offset=2.5)
    ax.set_xlabel("Residue number")
    ax.set_ylabel(r"RMSF ($\AA$)")
    #ax.set_ylim(top=9.5)

    #ax.axvspan(4, 9, alpha=0.1, color='red') # mark the hydophobic patch region
    # ax.axvspan(127, 132, alpha=0.5, color='red')
    # ax.axvspan(250, 255, alpha=0.5, color='red')

    plt.legend()
    plt.savefig('rmsf.svg', format='svg')
    plt.show()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Calcualtes the RMSD of a trajectory with a PBC corrected (and centered) trajectory")
    parser.add_argument("-c", dest="gro_file", type=str, help='the coordinate file [.gro]')
    parser.add_argument("-f", dest="xtc_file", type=str, help='a corrected trajectory file, pbc artifacts removed and the protein centered')
    parser.add_argument("-sel", dest="selection", type=str, help='the name of the gorup you wish to use for the RMSD calculation, e.g. "name CA"')
    parser.add_argument("-nc", dest="no_chains", type=int, help='the number of chains in your protein, if not supplied will run RMSD for all chains')
    parser.add_argument("-plot", dest="plot_type", type=str, help='specify whether this is a RMSD plot of one run or multiple runs, keywords are: "single" and "multiple"')
    parser.add_argument("-avg", dest="avg_plt", type=str, help='simple y or n, this will plot an avergae of the runs if "yes"', default="n")
    options = parser.parse_args()


    u = load_syst(options.gro_file, options.xtc_file)

    resids, rmsf_data = run_rmsf(u)

    plot_rmsf(rmsf_data, resids)




    # if options.plot_type == "multiple":
    #
    #     print("Generating multiple RMSD plots...")
    #     fig = plot_multiple(plot_average=options.avg_plt)

    # elif options.plot_type == "single":
    #
    #     # if the number of chains is supplied, calc RMSD for each
    #     if options.no_chains:
    #
    #         data, title = rmsd_calc_chains(coord_file=options.gro_file,
    #         traj_file=options.xtc_file, selection_phrase=options.selection,
    #         no_chains=options.no_chains)
    #
    #         print("Generating RMSD plot for each chain in system...")
    #         fig = plot_chains(data, title)
    #
    #     else:
    #
    #         data, title = rmsd_calc(coord_file=options.gro_file,
    #         traj_file=options.xtc_file, selection_phrase=options.selection)
    #
    #         print("Generating a single RMSD plot...")
    #         fig = plot(data=data, plot_name=title)
    #
    # else:
    #
    #     print("You need to specify a plot type!")
    #     print("Use python calc_rmsd.py -h for more info")