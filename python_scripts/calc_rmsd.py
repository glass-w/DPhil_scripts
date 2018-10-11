import MDAnalysis as mda
import MDAnalysis.analysis.rms as rms
import numpy as np
import pandas as pd
from glob import glob
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

print(" ")
print("###############")
print(" ")
print("Using the following versions: ")
print(" ")
print("MDAnalysis " + str(mda.__version__))
print("NumPy " + str(np.__version__))
print(" ")
print("###############")
print(" ")

def rmsd_calc(coord_file, traj_file, selection_phrase):
    ''' A function that calculates the RMSD of a user defined region
    and runs through a supplied trajectory
    '''

    u = mda.Universe(coord_file, traj_file)

    print ("Calculating RMSD...")

    rmsd = []

    u.trajectory[0]
    ref = u.select_atoms(selection_phrase).positions # reference frame

    for i in range(len(u.trajectory)):

        u.trajectory[i]

        bb = u.select_atoms(selection_phrase).positions # updated frames

        rmsd.append((u.trajectory.time, rms.rmsd(ref, bb)))

        if i == len(u.trajectory) / 4:
            print ("...25 %")
        elif i == len(u.trajectory) / 2:
            print ("...50 %")
        elif i == (3 * len(u.trajectory) / 4):
            print ("...75 %")
        elif i == len(u.trajectory):
            print ("...100 %")
            print ("Calculation done.")

    rmsd = np.array(rmsd)

    syst_name = os.path.split(coord_file)

    np.savetxt(syst_name[-1].split(".g")[0] + "_rmsd_data.csv", rmsd, delimiter=",")

    return rmsd, syst_name[-1]


def plot(data, plot_name):
    '''A function to plot calculated RMSD data and save the image according to the .gro file input'''

    ax = plt.subplot(111)
    #ax.set_title(r"RMSD of $\beta$3 Trimer")# - (system: " + plot_name + ")")
    #ax.set_title(r"RMSD of Nav1.7 $\alpha$ Subunit Model")
    ax.plot((data[:, 0]/1000), (data[:, 1]), 'blue', lw=1.5, label=r"$R_G$", alpha=0.9)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD ($\AA$)")
    ax.figure.savefig(plot_name.split(".g")[0] + "_RMSD.svg", format='svg')
    plt.draw()


def rmsd_calc_chains(coord_file, traj_file, selection_phrase, no_chains):

    ''' A function that calculates the RMSD for each chain of
    a user defined region and runs through a supplied trajectory
    '''

    u = mda.Universe(coord_file, traj_file)

    chain_labels = {'0':'A', '1':'B', '2':'C', '3':'D', '4':'E', '5':'F', '6':'G'}
    chain_list = ['Chain ' + chain_labels[str(i)] for i in range(no_chains)]

    rmsd_chains_dict = {key:[] for key in chain_list}

    print("Calculating RMSD...")

    chain_length = int(len(u.select_atoms("name CA")) / no_chains)

    print("chain length: ", str(chain_length))
    prot_pos = 0

    for i, chain in enumerate(rmsd_chains_dict):
        print("Chain = ", i)

        u.trajectory[0]
        ref = u.select_atoms(selection_phrase).positions[prot_pos:prot_pos + chain_length] # reference frame

        for j in range(len(u.trajectory)): # change to j, ts in enumerate ...

            u.trajectory[j]

            bb = u.select_atoms(selection_phrase).positions[prot_pos:prot_pos + chain_length] # updated frames

            rmsd_chains_dict['Chain ' + chain_labels[str(i)]].append((u.trajectory.time, rms.rmsd(ref, bb)))

        prot_pos += chain_length

    syst_name = os.path.split(coord_file)

    #np.savetxt(syst_name[-1].split(".g")[0] + "_rmsd_data.csv", rmsd, delimiter=",")

    return rmsd_chains_dict, syst_name[-1]

def plot_chains(data_dict, plot_name):

    '''A function to plot calculated RMSD data for multiple chains
    and save the image according to the .gro file input.

    If more than one chain is detected the plot will show an RMSD
    for each chain in the system.
    '''

    num_of_chains = len(data_dict)

    if num_of_chains > 1:

        fig, ax = plt.subplots(len(data_dict), sharex=True, sharey=True)

        col = ['#a6cee3','#1f78b4','#b2df8a']

        for i, chain in enumerate(data_dict):

            df = np.array(data_dict[chain])

            ax[i].plot((df[:, 0]/1000), (df[:, 1]), c=col[i], lw=1.5, label=str(chain), alpha=0.7)
            ax[i].legend()
            ax[i].set_ylim(ymax=10)

        fig.text(0.5, 0.04, 'Time (ns)', ha='center')
        fig.text(0.04, 0.5, "RMSD (" + r'$\AA$' + ')', va='center', rotation='vertical')

    # If only one chain in the system.
    else:

        fig, ax = plt.subplots()

        col = ['#a6cee3', '#1f78b4', '#b2df8a']

        df = np.array(data_dict['Chain A'])

        ax.plot((df[:, 0] / 1000), (df[:, 1]), c=col[0], lw=1.5, alpha=0.7)
        #ax.set_ylim(ymax=16)

        fig.text(0.5, 0.04, 'Time (ns)', ha='center')
        fig.text(0.04, 0.5, "RMSD (" + r'$\AA$' + ')', va='center', rotation='vertical')


    plt.savefig(plot_name.split(".g")[0] + "_RMSD.svg", format='svg')
    plt.show()


def plot_multiple(plot_average):

    files = glob("*.csv")

    ax = plt.subplot(111)
    #ax.set_title(r"RMSD of $\beta$3 Trimer")
    #ax.set_title(r"RMSD of Nav1.7 $\alpha$ Subunit Model")
    # Temporary fix!
    #mean = np.mean(np.array([split_all[0], split_all[1]]), axis=0)

    ### Plotting ###

    y_list = []
    avg = []

    # plot repeats

    print(files)
    from numpy import genfromtxt

    for i, repeat in enumerate(range(len(files))):

        plot_data = genfromtxt(files[repeat], delimiter=',')

        #print plot_data

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

        ax.plot((time / 1000), rmsd, lw=1, alpha=0.75, label="Run " + str(i))

        y_list.append(rmsd)

    # plot average line

    if plot_average == "y":

        for row in range(len(rmsd)):

            # take the same particle from each repeat
            particle_list = [item[row] for item in y_list]

            # take the average over each run
            avg.append(np.average(particle_list))

        #rolling_mean = running_mean(avg, 100)

        #ax.plot(time / 1000, avg[:], lw=0.5, alpha=1, color='black', label='Average')

        #ax.plot(time / 1000, avg, lw=0.5, alpha=1, color='black', label='Average')

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD ($\AA$)")
    ax.legend()
    plt.title("RMSD")
    plt.draw()

    ax.figure.savefig("multiple_RMSD_line_plot.svg", format='svg', dpi=300)

    print("Plot complete.")

    plt.clf()

    ax = plt.subplot(111)
    y_list = []

    ############################################
    # plot histogram repeats

    for i, repeat in enumerate(range(len(files))):

        plot_data = genfromtxt(files[repeat], delimiter=',')

        # print plot_data

        time = plot_data[:, 0]
        rmsd = plot_data[:, 1]

        sns.kdeplot(rmsd, label="Run " + str(i))

    ax.set_xlabel("RMSD ($\AA$)")
    ax.set_ylabel("Count")
    ax.legend()
    plt.title("Histogram of RMSD Values")
    plt.draw()

    ax.figure.savefig("multiple_RMSD_histogram_plot.svg", format='svg', dpi=300)

    #############################################
    # plot violin plot of repeats
    plt.clf()

    from collections import OrderedDict

    df = OrderedDict()

    #df = {}

    for i, repeat in enumerate(range(len(files))):

        run_i = 'Run {}'.format(i)

        plot_data = genfromtxt(files[repeat], delimiter=',')


        time = plot_data[:, 0]
        rmsd = plot_data[:, 1]

        df[run_i] = rmsd

        #sns.kdeplot(rmsd, label="Run " + str(i))

    print(df)

    df = pd.DataFrame.from_dict(df)

    print(df)
    #
    sns.violinplot(data=df)
    #ax.set_xlabel("RMSD ($\AA$)")
    ax.set_ylabel("RMSD ($\AA$)")
    # ax.legend()
    # plt.title("Histogram of RMSD Values")
    plt.show()
    #
    # ax.figure.savefig("multiple_RMSD_histogram_plot.svg", format='svg', dpi=300)


def running_mean(x, N):

    cumsum = np.cumsum(np.insert(x, 0, 0))

    return (cumsum[N:] - cumsum[:-N]) / N


# To Do: add an average line to the multiple plots...

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Calcualtes the RMSD of a trajectory with a PBC corrected (and centered) trajectory")
    parser.add_argument("-c", dest="gro_file", type=str, help='the coordinate file [.gro]')
    parser.add_argument("-f", dest="xtc_file", type=str, help='a corrected trajectory file, pbc artifacts removed and the protein centered')
    parser.add_argument("-sel", dest="selection", type=str, help='the name of the gorup you wish to use for the RMSD calculation, e.g. "name CA"')
    parser.add_argument("-nc", dest="no_chains", type=int, help='the number of chains in your protein, if not supplied will run RMSD for all chains')
    parser.add_argument("-plot", dest="plot_type", type=str, help='specify whether this is a RMSD plot of one run or multiple runs, keywords are: "single" and "multiple"')
    parser.add_argument("-avg", dest="avg_plt", type=str, help='simple y or n, this will plot an avergae of the runs if "yes"', default="n")
    options = parser.parse_args()

    if options.plot_type == "multiple":

        print("Generating multiple RMSD plots...")
        fig = plot_multiple(plot_average=options.avg_plt)

    elif options.plot_type == "single":

        # if the number of chains is supplied, calc RMSD for each
        if options.no_chains:

            data, title = rmsd_calc_chains(coord_file=options.gro_file,
            traj_file=options.xtc_file, selection_phrase=options.selection,
            no_chains=options.no_chains)

            print("Generating RMSD plot for each chain in system...")
            fig = plot_chains(data, title)

        else:

            data, title = rmsd_calc(coord_file=options.gro_file,
            traj_file=options.xtc_file, selection_phrase=options.selection)

            print("Generating a single RMSD plot...")
            fig = plot(data=data, plot_name=title)

    else:

        print("You need to specify a plot type!")
        print("Use python calc_rmsd.py -h for more info")
