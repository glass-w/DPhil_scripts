import MDAnalysis as mda
import MDAnalysis.analysis.rms as rms
import numpy as np
from glob import glob
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

print " "
print "###############"
print " "
print "Using the following versions: "
print " "
print "MDAnalysis " + str(mda.__version__)
print "NumPy " + str(np.__version__)
print "Pandas " + str(pd.__version__)
print " "
print "###############"
print " "

def rmsd_calc(coord_file, traj_file, selection_phrase):
    ''' A function that calculates the RMSD of a user defined region
    and runs through a supplied trajectory
    '''

    u = mda.Universe(coord_file, traj_file)

    print ("Calculating RMSD...")

    rmsd = []

    u.trajectory[1]
    ref = u.select_atoms(selection_phrase).positions # reference frame

    for i in range(len(u.trajectory)):

        u.trajectory[i]

        bb = u.select_atoms(selection_phrase).positions # updated frames

        rmsd.append((u.trajectory.time, rms.rmsd(ref, bb)))

        i += 1

        if i == len(u.trajectory)/4:
            print ("...25 %")
        elif i == len(u.trajectory)/2:
            print ("...50 %")
        elif i == len(u.trajectory)*(3/4):
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
    ax.set_title(r"RMSD of $\beta$3 Trimer EC Domain - (system: " + plot_name + ")")
    ax.plot((data[:, 0]/1000), (data[:, 1]), 'blue', lw=1.5, label=r"$R_G$", alpha=1)
    ax.set_xlabel("Time / ns")
    ax.set_ylabel("RMSD from t = 0 / $\AA$")
    ax.figure.savefig(plot_name.split(".g")[0] + "_RMSD.svg", format='svg')
    plt.draw()

def plot_multiple():

    files = glob("./*.csv")

    print files

    dataLists =[]

    for file in files:

        dataLists[file.split(".")[0]].append(pd.read_csv(file))


    ax = plt.subplot(111)
    ax.set_title(r"RMSD of $\beta$3 Trimer EC Domain")

    # Temporary fix!
    mean = np.mean(np.array([split_all[0], split_all[1]]), axis=0)

    ### Plotting ###

    for repeat in range(len(files)):

        plot_data = np.array(pd.DataFrame(pd.read_csv(files[repeat])))

        ax.plot((plot_data[:, 0] / 1000), (plot_data[:, 1]), lw=1, alpha=0.25,
                label='Run ' + str(repeat))

    ax.plot(mean[:, 0] / 1000, mean[:, 1], lw=0.5, alpha=1, color='black', label='Average')

    ax.set_xlabel("Time / ns")
    ax.set_ylabel("RMSD from t = 0 / $\AA$")
    ax.legend()
    plt.draw()

    ax.figure.savefig("multiple_RMSD.svg", format='svg')

    print "Plot complete."


# To Do: add an average line to the multiple plots...

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Calcualtes the RMSD of a trajectory with a PBC corrected (and centered) trajectory")
    parser.add_argument("-c", dest="gro_file", type=str, help='the coordinate file [.gro]')
    parser.add_argument("-f", dest="xtc_file", type=str, help='a corrected trajectory file, pbc artifacts removed and the protein centered')
    parser.add_argument("-sel", dest="selection", type=str, help='the name of the gorup you wish to use for the RMSD calculation, e.g. "name CA"')
    parser.add_argument("-plot", dest="plot_type", type=str, help='specify whether this is a RMSD plot of one run or multiple runs, keywords are: "single" and "multiple"')
    options = parser.parse_args()

    if options.plot_type == "multiple":

        print "Generating multiple RMSD plots..."
        fig = plot_multiple()

    elif options.plot_type == "single":

        data, title = rmsd_calc(coord_file=options.gro_file, traj_file=options.xtc_file, selection_phrase=options.selection)
        fig = plot(data=data, plot_name=title)

    else:

        print "You need to specify a plot type!"
        print "Use python calc_rmsd.py -h for more info"




