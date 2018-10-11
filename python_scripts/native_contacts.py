import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(color_codes=True)
import sys
import argparse


def load_uni(gro_file, xtc_file):

    u = mda.Universe(gro_file, xtc_file)

    return u

def plot_single(data):

    average_contacts = np.mean(data)

    f, ax = plt.subplots()


    # Histogram

    ax.hist(ca.timeseries[:, 1], bins=50)

    ax.set(xlabel='fraction of native contacts', ylabel='count')

    plt.savefig('salt_bridge_native_contacts_hist.svg', format='svg', dpi=300)
    plt.show()

    plt.clf()

    # Line Plot

    ax = sns.lineplot(ca.timeseries[:, 0], ca.timeseries[:, 1])

    ax.set(xlabel='Frame', ylabel='fraction of native contacts')

    plt.savefig('salt_bridge_native_contacts_line.svg', format='svg', dpi=300)

    plt.show()


def plot_multiple(data_list):

    plot_data = []

    # load the numpy arrays in
    for array in data_list:

        plot_data.append(np.load(array))


    f, ax = plt.subplots()

    for i, array in enumerate(plot_data):

        ax = sns.kdeplot(array, shade=True, label="Run " + str(i))

    ax.set(xlabel='fraction of native contacts', ylabel='count')
    plt.legend()

    plt.savefig('salt_bridge_contacts_multiple_systems_hist.svg', format='svg', dpi=300)
    plt.show()

    plt.clf()

    for i, array in enumerate(plot_data):

        ax = sns.lineplot(array, label="Run " + str(i))

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculates the native contacts")
    parser.add_argument("-c", dest="gro_file", help='the coordinate file [.gro]')
    parser.add_argument("-f", dest="xtc_file", help='a corrected trajectory file, pbc artifacts removed and the protein centered')

    parser.add_argument("-dl", dest="data_list", type=str, nargs='+', help='a list of numpy arrays (.npy) that have been calculated previously.')

    parser.add_argument("-pt", dest="plot_type", type=str, help="type of plot, either 'single' or 'multiple'")

    parser.add_argument("-cut_off", dest="cut_off", type=float, help="the cut-off value")

    options = parser.parse_args()


    if options.plot_type == "single":

        uni = load_uni(options.gro_file, options.xtc_file)

        sel_basic = "(resname ARG LYS) and (name NH* NZ)"
        sel_acidic = "(resname ASP GLU) and (name OE* OD*)"

        acidic = uni.select_atoms(sel_acidic)
        basic = uni.select_atoms(sel_basic)

        ca = contacts.Contacts(uni, selection=(sel_acidic, sel_basic), refgroup=(acidic, basic), radius=options.cut_off,
                               verbose=True)

        ca.run()

        np.save('contact_analysis_data.npy', ca.timeseries[:, 1])

        plot_single(ca.timeseries[:, 1])

    elif options.plot_type == "multiple":

        plot_multiple(options.data_list)

    else:

        print("You have not specified a plot type!")



