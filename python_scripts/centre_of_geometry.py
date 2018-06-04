import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import numpy as np
from scipy.interpolate import spline
from matplotlib import pyplot as plt
import argparse
import sys


def area(a, b, c):
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a))


def cog(gro_file, traj_file):
    print gro_file
    print traj_file

    u = mda.Universe(gro_file, traj_file)

    # define regions of which you want to calculate the centre of geometry
    # ig_a = u.select_atoms("resid 1:123")
    # ig_b = u.select_atoms("resid 167:290")
    # ig_c = u.select_atoms("resid 333:457")

    # ig_a = u.select_atoms("resid 1:123")
    # ig_b = u.select_atoms("resid 124:246")
    # ig_c = u.select_atoms("resid 247:369")

    ig_a = u.select_atoms("resid 13:123")
    ig_b = u.select_atoms("resid 136:246")
    ig_c = u.select_atoms("resid 259:369")

    r1_store, r2_store, r3_store = [], [], []

    d1_store, d2_store, d3_store = [], [], []

    area_store = []
    time = []

    for frame in u.trajectory[::10]:

        time.append((u.trajectory.time) / 1000)
        print "Time = " + str(((u.trajectory.time) / 1000)) + " ns"

        ig_a_com = ig_a.select_atoms("name CA").center_of_mass()
        ig_b_com = ig_b.select_atoms("name CA").center_of_mass()
        ig_c_com = ig_c.select_atoms("name CA").center_of_mass()

        cent_of_prots = (ig_a.select_atoms("name CA") + ig_b.select_atoms("name CA") + ig_c.select_atoms("name CA")).center_of_mass()

        #Get distances between each protein and the CoM of all of them

        d1 = np.linalg.norm(ig_a_com - cent_of_prots)
        d1_store.append(d1)
        d2 = np.linalg.norm(ig_b_com - cent_of_prots)
        d2_store.append(d2)
        d3 = np.linalg.norm(ig_c_com - cent_of_prots)
        d3_store.append(d3)

        # Get distances between CoMs
        r1 = np.linalg.norm(ig_a_com - ig_c_com)
        r1_store.append(r1)

        r2 = np.linalg.norm(ig_a_com - ig_b_com)
        r2_store.append(r2)

        r3 = np.linalg.norm(ig_b_com - ig_c_com)
        r3_store.append(r3)

        # Get vectors between CoMs

        vec1 = np.array([ig_a_com[0], ig_a_com[1]])
        vec2 = np.array([ig_b_com[0], ig_b_com[1]])
        vec3 = np.array([ig_c_com[0], ig_c_com[1]])

        tria_area = area(vec1, vec2, vec3)

        area_store.append(tria_area)

    ##################
    #                #
    #    PLOTTING    #
    #                #
    ##################

    xnew = np.linspace(np.min(time), np.max(time), 150)  # No. represents number of points to make between min and max

    c = ['#a6cee3', '#1f78b4', '#b2df8a']
    c2 = ['#1b9e77','#d95f02','#7570b3']

    ##### Plot distances between CoM of each protein and the CoM of the protein cluster #####

    power_smooth1 = spline(time, d1_store, xnew)
    power_smooth2 = spline(time, d2_store, xnew)
    power_smooth3 = spline(time, d3_store, xnew)

    ax = plt.subplot(111)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Distance from Cluster Centre ($\AA$)")

    ax.plot(xnew, power_smooth1, c2[0], label='Chain A', alpha=0.9)
    ax.plot(time, d1_store, c2[0], alpha=0.2)

    ax.plot(xnew, power_smooth2, c2[1], label='Chain B', alpha=0.9)
    ax.plot(time, d2_store, c2[1], alpha=0.2)

    ax.plot(xnew, power_smooth3, c2[2], label='Chain C', alpha=0.9)
    ax.plot(time, d3_store, c2[2], alpha=0.2)

    if options.plot_avg_line_flag == True:

        ax.plot(xnew, np.mean([power_smooth1, power_smooth2, power_smooth3], axis=0), c='black', label='Average', alpha=1)

    ax.legend()
    plt.title("Distance of " + r'$\beta$3' + "Trimer EC Domains & Trimer Centre (C4-C26 Broken)")
    plt.savefig(str(options.file_name) + '_distance_from_cent.svg', dpi=300)
    plt.clf()


    ##### Plot distances between CoM's #####

    power_smooth1 = spline(time, r1_store, xnew)
    power_smooth2 = spline(time, r2_store, xnew)
    power_smooth3 = spline(time, r3_store, xnew)

    ax = plt.subplot(111)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Distance ($\AA$)")

    ax.plot(xnew, power_smooth1, c[0], label=('d' + r'$_{A - C}$'), alpha=0.9)
    ax.plot(time, r1_store, c[0], alpha=0.2)

    ax.plot(xnew, power_smooth2, c[1], label=('d' + r'$_{A - B}$'), alpha=0.9)
    ax.plot(time, r2_store, c[1], alpha=0.2)

    ax.plot(xnew, power_smooth3, c[2], label=('d' + r'$_{B - C}$'), alpha=0.9)
    ax.plot(time, r3_store, c[2], alpha=0.2)

    if options.plot_avg_line_flag == True:

        ax.plot(xnew, np.mean([power_smooth1, power_smooth2, power_smooth3], axis=0), c='black', label='Average', alpha=1)

    ax.legend()
    plt.title("Distance Between EC Domains of " + r'$\beta$3' + " (water box, C4-C26 Broken)")
    plt.savefig(str(options.file_name) + '_distance.svg', dpi=300)
    plt.clf()

    #######################################

    #### Plot triangle area ####

    ax = plt.subplot(111)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Area ($\AA ^2$)")
    plt.title("Area Between EC Domains of " + r'$\beta$3' + " (water box, C4-C26 Broken)")

    ax.plot(xnew, spline(time, area_store, xnew))

    plt.savefig(str(options.file_name) + '_area.svg', dpi=300)

    plt.clf()

    ###########################

    ##################

    # Write out arrays

    dist_cent_data = np.array([np.array(d1_store), np.array(d2_store), np.array(d3_store)])

    distance_data = np.array([np.array(r1_store), np.array(r2_store), np.array(r3_store)])
    area_data = np.array([area_store])
    time_data = np.array([time])

    np.save(str(options.file_name) + '_dist_cent_data', dist_cent_data)
    np.save(str(options.file_name) + '_distance_data', distance_data)
    np.save(str(options.file_name) + '_area_data', area_data)
    np.save(str(options.file_name) + '_time_data', time)

    #################

    return distance_data, area_data, time_data

def plot_avg(time, *args):

    time = np.load(time)

    c1 = ['#49006a', '#ae017e', '#f768a1', '#fcc5c0']
    #c1 =['#800026', '#e31a1c', '#fd8d3c', '#fed976']

    xnew = np.linspace(np.min(time), np.max(time), 150)

    print args[0]
    area_data_dict = {}

    i = 0

    for arg in args[0]:

        area_data_dict["df{0}".format(i)] = np.load(str(arg))

        i += 1

    ax = plt.subplot(111)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Area ($\AA ^2$)")
    plt.title("Area Enclosed by EC Domains of " + r'$\beta$3' + " (water box, CYS4-26 broken)")
    ax.set_ylim([260, 380])

    j = 0

    for df in area_data_dict:

        ax.plot(xnew, spline(time, area_data_dict[df][0], xnew), alpha=0.4, color=c1[j], label='Run ' + str(j))

        j += 1

    area_mean = np.mean([area_data_dict["df0"][0], area_data_dict["df1"][0], area_data_dict["df2"][0]], axis=0)#, area_data_dict["df3"][0]], axis=0)

    #print area_mean

    ax.plot(xnew, spline(time, area_mean, xnew), alpha=1, color='black', label='Average')

    ax.legend()

    if len(options.file_name) > 0:

        plt.savefig(str(options.file_name) + '_area.svg', dpi=300)
        plt.clf()

    else:

        plt.savefig('area_plots.svg', dpi=300)
        plt.clf()

def native_contacts(gro_list, traj_list, num_proteins, N_resid):
    c_pal = ['#225ea8', '#41b6c4', '#c7e9b4']

    ax = plt.subplot(111)
    ax.set(xlabel='Frame', ylabel='% of Native Contacts', title='Native Contacts')
    ax.set_ylim([0, 100])

    if len(gro_list) == 1:
        compare = False

        u = mda.Universe(gro_list[0], traj_list[0])
        protein_length = len(u.select_atoms("protein and name CA"))

        print "One universe loaded..."

        # Initialise dicts, first for each protein, second for the residues you wish to calc the Native Conts for.

        protein_dict = {}
        protein_res_to_calc_dict = {}
        regions_for_NC_analysis = {}

        # Get protein information and store in dictionaries

        for j in range(num_proteins):

            print "j = ", j

            if j == 0:

                start = ((j * (protein_length / num_proteins)) + 1)

                end = (j * (protein_length / num_proteins)) + N_resid

                print "start and end = ", start, end

                protein_dict["p" + str(j)] = u.select_atoms("resid " + str(1) + ":" + str(N_resid))

                protein_res_to_calc_dict["p" + str(j) + "res"] = {}

                protein_res_to_calc_dict["p" + str(j) + "res"]["sel"] = u.select_atoms("(resid " + str(1) + ":" + str(N_resid)
                                                                                + ") and (not name CA N O H*)")


                print protein_res_to_calc_dict["p" + str(j) + "res"]["sel"].resnames
                print protein_res_to_calc_dict["p" + str(j) + "res"]["sel"].resids

                # print protein_res_to_calc_dict
                # store the string used to select this part of the protein for future use

                protein_res_to_calc_dict["p" + str(j) + "res"]["sel_string"] = "resid " + str(1) + ":" + str(N_resid)

                j += 1

            elif j >= 1:

                # make sure you're in the right place for the next protein

                # print "one protein length", protein_length / num_proteins
                # print "num proteins", num_proteins
                # print "n resid", N_resid

                start = (j * (protein_length / num_proteins)) + 1

                end = (j * (protein_length / num_proteins)) + N_resid

                print "protein length = ", protein_length
                print "num prots = ", num_proteins

                print "start and end = ", start, end

                protein_res_to_calc_dict["p" + str(j) + "res"] = {}

                # populate dicts
                protein_dict["p" + str(j)] = u.select_atoms("resid " + str(start) + ":" +
                                                            str(end))

                protein_res_to_calc_dict["p" + str(j) + "res"]["sel"] = u.select_atoms("(resid " + str(start) + ":"
                                                                                + str(end)
                                                                                + ") and (not name CA N O H*)")

                print protein_res_to_calc_dict["p" + str(j) + "res"]["sel"].resnames
                print protein_res_to_calc_dict["p" + str(j) + "res"]["sel"].resids

                # store the string used to select this part of the protein for future use


                protein_res_to_calc_dict["p" + str(j) + "res"]["sel_string"] = "resid " + str(start) + ":" + str(end)

                j += 1

        print protein_dict

        for k, key in enumerate(protein_dict):

            print "current key = ", key

            # A
            sel0 = protein_res_to_calc_dict["p0res"]["sel_string"]  # + " and (not name CA N O H*)"
            # sel0 = "(resid " + str(1) + ":" + str(N_resid) + ") and (not name CA N O H*)"
            ref0 = protein_res_to_calc_dict["p0res"]["sel"]

            # B
            sel1 = protein_res_to_calc_dict["p1res"]["sel_string"]
            ref1 = protein_res_to_calc_dict["p1res"]["sel"]

            # C
            sel2 = protein_res_to_calc_dict["p2res"]["sel_string"]
            ref2 = protein_res_to_calc_dict["p2res"]["sel"]

            # For first protein
            if key == "p0":

                comb_ref = ref1 + ref2

                print ""
                print "selections"
                print ""
                print "sel0 = ", sel0
                print "sel1 = ", sel1
                print "sel2 = ", sel2
                print ""

                print "references"
                print ""
                print "ref0 = ", len(ref0)
                print "ref1 = ", len(ref1)
                print "ref2 = ", len(ref2)
                print ""

                print "selection check"
                print sel0 + " and not name CA N O H*"
                print sel1 + " or " + sel2 + " and not name CA N O H*"

                #"(resid 1:10) and (not name CA N O H*)"
                #"(resid 124:134 or resid 247:257) and (not name CA N O H*)"

                # A - BC
                regions_for_NC_analysis["A_BC"] = contacts.Contacts(u, selection=(sel0 + " and not name CA N O H*", sel1 + " or " + sel2 + " and not name CA N O H*"),
                                                                          refgroup=(ref0, comb_ref), radius=6.0)

                regions_for_NC_analysis["A_BC"].run()

                print np.mean(regions_for_NC_analysis["A_BC"].timeseries[:, 1])

                plt.plot(regions_for_NC_analysis["A_BC"].timeseries[:, 0],
                         (100 * regions_for_NC_analysis["A_BC"].timeseries[:, 1]), alpha=0.5, label='Chain A', c=c_pal[k])

            elif key == "p1":

                comb_ref = ref0 + ref2

                # B - AC
                regions_for_NC_analysis["B_AC"] = contacts.Contacts(u, selection=(sel1 + " and not name CA N O H*", sel0 + " or " + sel2 + " and not name CA N O H*"),
                                                                    refgroup=(ref1, comb_ref), radius=6.0)

                regions_for_NC_analysis["B_AC"].run()

                print np.mean(regions_for_NC_analysis["B_AC"].timeseries[:, 1])

                plt.plot(regions_for_NC_analysis["B_AC"].timeseries[:, 0],
                         (100 * regions_for_NC_analysis["B_AC"].timeseries[:, 1]), alpha=0.5, label='Chain B', c=c_pal[k])

            elif key =="p2":

                comb_ref = ref0 + ref1

                # C - AB
                regions_for_NC_analysis["C_AB"] = contacts.Contacts(u, selection=(sel2 + " and not name CA N O H*", sel0 + " or " + sel1 + " and not name CA N O H*"),
                                                                    refgroup=(ref2, comb_ref), radius=6.0)

                regions_for_NC_analysis["C_AB"].run()

                print np.mean(regions_for_NC_analysis["C_AB"].timeseries[:, 1])

                plt.plot(regions_for_NC_analysis["C_AB"].timeseries[:, 0],
                         (100 * regions_for_NC_analysis["C_AB"].timeseries[:, 1]), alpha=0.5, label='Chain C', c=c_pal[k])

            else:
                continue

        # for key in regions_for_NC_analysis:

        avg = 100 * (np.mean([regions_for_NC_analysis["A_BC"].timeseries[:, 1], regions_for_NC_analysis["B_AC"].timeseries[:, 1],
                       regions_for_NC_analysis["C_AB"].timeseries[:, 1]], axis=0))

        #print ({k: np.mean([v[:, 1]], axis=0) for k, v in regions_for_NC_analysis.items()})

        print np.mean([regions_for_NC_analysis[k].timeseries[:, 1] for k, v in regions_for_NC_analysis.items()], axis=0)

        plt.plot(avg, color='black', alpha=1, label='Avg')

        plt.legend()

        plt.savefig('NC_plot_' + str(num_proteins) + "prot_syst" + '.svg', format='svg')

        plt.show()
        plt.clf()

        ###### COMPARE TWO SYSTEMS #######

    elif len(gro_list) != 1:
        compare = True

        u_list = [mda.Universe(gro_list[0], traj_list[0]), mda.Universe(gro_list[1], traj_list[1])]

        print "Two universes loaded..."

        # Initialise dicts, first for each protein, second for the residues you wish to calc the Native Conts for.

        protein_dict = {}
        protein_res_to_calc_dict = {}
        regions_for_NC_analysis = {}
        avg_dict = {}

    # Get protein information and store in dictionaries, for each system

        # Initialise dicts for each system present

        for s, system in enumerate(gro_list):

            protein_dict["s" + str(s)] = {}
            protein_res_to_calc_dict["s" + str(s)] = {}
            regions_for_NC_analysis["s" + str(s)] = {}

        for s, system in enumerate(gro_list):

            print s

            if s > 0:
                c_pal = ['#ae017e','#f768a1', '#fcc5c0']

            else:
                c_pal = ['#225ea8', '#41b6c4', '#c7e9b4']

            for j in range(num_proteins):

                if j == 0:

                    #system, protein, region ----> dict structure

                    protein_dict["s" + str(s)]["p" + str(j)] = u_list[s].select_atoms("resid " + str(1) + ":" + str(N_resid))

                    protein_res_to_calc_dict["s" + str(s)]["p" + str(j) + "res"] = {}

                    protein_res_to_calc_dict["s" + str(s)]["p" + str(j) + "res"]["sel"] = u_list[s].select_atoms(
                        "(resid " + str(1) + ":" + str(N_resid)
                        + ") and (not name CA N O H*)")

                    # print protein_res_to_calc_dict["s" + str(s)]["p" + str(j) + "res"]["sel"].resnames
                    # print protein_res_to_calc_dict["s" + str(s)]["p" + str(j) + "res"]["sel"].resids

                    # print protein_res_to_calc_dict
                    # store the string used to select this part of the protein for future use

                    protein_res_to_calc_dict["s" + str(s)]["p" + str(j) + "res"]["sel_string"] = "resid " + str(1) + ":" + str(N_resid)

                    j += 1

                elif j >= 1:

                    # make sure you're in the right place for the next protein, for the current system

                    protein_length = len(u_list[s].select_atoms("protein and name CA"))

                    start = (j * (protein_length / num_proteins)) + 1

                    end = (j * (protein_length / num_proteins)) + N_resid

                    protein_res_to_calc_dict["s" + str(s)]["p" + str(j) + "res"] = {}

                    # populate dicts
                    protein_dict["s" + str(s)]["p" + str(j)] = u_list[s].select_atoms("resid " + str(start) + ":" +
                                                                str(end))

                    protein_res_to_calc_dict["s" + str(s)]["p" + str(j) + "res"]["sel"] = u_list[s].select_atoms("(resid " + str(start) + ":"
                                                                                           + str(end)
                                                                                           + ") and (not name CA N O H*)")

                    # print protein_res_to_calc_dict["s" + str(s)]["p" + str(j) + "res"]["sel"].resnames
                    # print protein_res_to_calc_dict["s" + str(s)]["p" + str(j) + "res"]["sel"].resids

                    # store the string used to select this part of the protein for future use

                    protein_res_to_calc_dict["s" + str(s)]["p" + str(j) + "res"]["sel_string"] = "resid " + str(start) + ":" + str(end)

                    j += 1

            for k, key in enumerate(protein_dict["s" + str(s)]):

                print "current system, chain = ", s, key

                # A
                sel0 = protein_res_to_calc_dict["s" + str(s)]["p0res"]["sel_string"]  # + " and (not name CA N O H*)"
                # sel0 = "(resid " + str(1) + ":" + str(N_resid) + ") and (not name CA N O H*)"
                ref0 = protein_res_to_calc_dict["s" + str(s)]["p0res"]["sel"]

                # B
                sel1 = protein_res_to_calc_dict["s" + str(s)]["p1res"]["sel_string"]
                ref1 = protein_res_to_calc_dict["s" + str(s)]["p1res"]["sel"]

                # C
                sel2 = protein_res_to_calc_dict["s" + str(s)]["p2res"]["sel_string"]
                ref2 = protein_res_to_calc_dict["s" + str(s)]["p2res"]["sel"]

                # For first protein
                if key == "p0":

                    comb_ref = ref1 + ref2

                    # "(resid 1:10) and (not name CA N O H*)"
                    # "(resid 124:134 or resid 247:257) and (not name CA N O H*)"

                    # A - BC
                    regions_for_NC_analysis["s" + str(s)]["A_BC"] = contacts.Contacts(u_list[s], selection=(
                    sel0 + " and not name CA N O H*", sel1 + " or " + sel2 + " and not name CA N O H*"),
                                                                        refgroup=(ref0, comb_ref), radius=6.0)

                    regions_for_NC_analysis["s" + str(s)]["A_BC"].run()

                    print np.mean(regions_for_NC_analysis["s" + str(s)]["A_BC"].timeseries[:, 1])

                    plt.plot(regions_for_NC_analysis["s" + str(s)]["A_BC"].timeseries[:, 0],
                             (100 * regions_for_NC_analysis["s" + str(s)]["A_BC"].timeseries[:, 1]), alpha=0.3, label='Chain A (system ' + str(s) + ')',
                             c=c_pal[k])

                elif key == "p1":

                    comb_ref = ref0 + ref2

                    # B - AC
                    regions_for_NC_analysis["s" + str(s)]["B_AC"] = contacts.Contacts(u_list[s], selection=(
                    sel1 + " and not name CA N O H*", sel0 + " or " + sel2 + " and not name CA N O H*"),
                                                                        refgroup=(ref1, comb_ref), radius=6.0)

                    regions_for_NC_analysis["s" + str(s)]["B_AC"].run()

                    print np.mean(regions_for_NC_analysis["s" + str(s)]["B_AC"].timeseries[:, 1])

                    plt.plot(regions_for_NC_analysis["s" + str(s)]["B_AC"].timeseries[:, 0],
                             (100 * regions_for_NC_analysis["s" + str(s)]["B_AC"].timeseries[:, 1]), alpha=0.3, label='Chain B (system ' + str(s) + ')',
                             c=c_pal[k])

                elif key == "p2":

                    comb_ref = ref0 + ref1

                    # C - AB
                    regions_for_NC_analysis["s" + str(s)]["C_AB"] = contacts.Contacts(u_list[s], selection=(
                    sel2 + " and not name CA N O H*", sel0 + " or " + sel1 + " and not name CA N O H*"),
                                                                        refgroup=(ref2, comb_ref), radius=6.0)

                    regions_for_NC_analysis["s" + str(s)]["C_AB"].run()

                    print np.mean(regions_for_NC_analysis["s" + str(s)]["C_AB"].timeseries[:, 1])

                    plt.plot(regions_for_NC_analysis["s" + str(s)]["C_AB"].timeseries[:, 0],
                             (100 * regions_for_NC_analysis["s" + str(s)]["C_AB"].timeseries[:, 1]), alpha=0.3, label='Chain C (system ' + str(s) + ')',
                             c=c_pal[k])

                else:
                    continue

            time_data = np.arange(0, 2501, 1) # need to generalise this

            avg = np.mean([regions_for_NC_analysis["s" + str(s)]["A_BC"].timeseries[:, 1],
                                  regions_for_NC_analysis["s" + str(s)]["B_AC"].timeseries[:, 1],
                                  regions_for_NC_analysis["s" + str(s)]["C_AB"].timeseries[:, 1]], axis=0)

            avg_spl = 100 * avg

            xnew = np.linspace(np.min(time_data), np.max(time_data), 300)  # No. represents number of points to make between min and max

            power_smooth_avg = spline(time_data, avg_spl, xnew)

            print time_data
            print avg
            print power_smooth_avg
            print ""

            print time_data[0]
            print avg[0]
            print power_smooth_avg[0]

            plt.plot(xnew, power_smooth_avg, label='Avg (system ' + str(s) + ')', alpha=1, color='black')

            #plt.plot(avg, color='black', alpha=1, label='Avg (system ' + str(s) + ')')

        # need to add this method below for a more general calc of mean -- ???
        #print np.mean([regions_for_NC_analysis[k].timeseries[:, 1] for k, v in regions_for_NC_analysis["s" + str(s)].items()], axis=0)

        plt.legend()

        plt.savefig('compare_plot_' + str(num_proteins) + "prot_syst" + '.svg', format='svg')

        plt.show()
        plt.clf()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculates the centre of geometry between three defined protein regions")
    parser.add_argument("-c", dest="gro_file", type=str, nargs='+', help='the coordinate file [.gro]')
    parser.add_argument("-f", dest="xtc_file", type=str, nargs='+',
                        help='a corrected trajectory file, pbc artifacts removed and the protein centered')
    parser.add_argument("-avg", dest="plot_avg_line_flag", action="store_true", default=False, help='include if you want an average line for the distance plots')

    parser.add_argument("-native_c", dest="native_c", action="store_true", help='to carry out Native contact analysis')

    parser.add_argument("-n", dest="file_name", type=str, help="the name of the output files")
    parser.add_argument("-plota", dest="plot_avg_flag", action="store_true", default=False, help='include if you want to plot avergaes of data already collected')
    parser.add_argument("-areafiles", dest="files_for_avg", nargs='+', help='a list of the files containing area info')
    parser.add_argument("-timefile", dest="time_file_for_avg", help='an array containing the time for plotting')

    options = parser.parse_args()

    if options.native_c == True:

        native_contacts(options.gro_file, options.xtc_file, num_proteins=3, N_resid=10)

    elif options.plot_avg_flag == False:

        cog(options.gro_file, options.xtc_file)

    elif options.plot_avg_flag == True:

        plot_avg(options.time_file_for_avg, options.files_for_avg)


