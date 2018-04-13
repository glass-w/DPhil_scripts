from __future__ import print_function
from __future__ import print_function
from __future__ import print_function
from sys import stdout
import time
import argparse
from random import randint
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as plt
import MDAnalysis as mda

print ' '
print "### Modules Used ####"
print ' '
print "NumPy Version: " + np.__version__
print "MDAnalysis Version: " + mda.__version__
print ' '
print "####################"
print ' '

'''

A script to analyse the diffusion of proteins within a lipid bilayer

Author: W. Glass (with reference to script by Tyler Reddy)

See: https://github.com/tylerjereddy/diffusion_analysis_MD_simulations/blob/master/README.md

'''


def centroid_array_production_protein(protein_sel, num_protein_copies):
    """ Creates a dictionary containing the "centre" of each protein in the system (i.e. a centroid)"""

    dictionary_centroid_arrays = {}

    full_protein_coord_array = protein_sel.positions

    list_individual_protein_coord_arrays = np.split(full_protein_coord_array, num_protein_copies)

    list_per_protein_centroids = [np.average(protein_coord_array, axis=0) for protein_coord_array in
                                  list_individual_protein_coord_arrays]

    dictionary_centroid_arrays['protein'] = np.array(list_per_protein_centroids)

    return dictionary_centroid_arrays  # Now have a dictionary of the 'centre' of each protein in the membrane


def msd(grofile, trajfile):

    def increments(u):

        num_of_elements = int(round(len(u.trajectory)) / 3)

        frst_qrtr = np.linspace(start=1, stop=num_of_elements / 4, num=num_of_elements / 4)
        scnd_qrtr = np.linspace(start=(num_of_elements / 4) + 1, stop=(num_of_elements / 2) - 1,
                                num=num_of_elements / 16)
        last_half = np.linspace(start=num_of_elements / 2, stop=num_of_elements, num=(num_of_elements / 32) + 1)

        print(np.concatenate((np.round(frst_qrtr), np.round(scnd_qrtr), np.round(last_half))))

        return np.concatenate((np.round(frst_qrtr), np.round(scnd_qrtr), np.round(last_half)))

    # def make_number_list(Universe):
    #
    #     # total_time = u.trajectory[-1].time
    #     #
    #     # time_counter = 0
    #     #
    #     # for frame in u.trajectory:
    #     #
    #     #     time_counter += 1
    #     #
    #     #     if round(u.trajectory.time) >= round(total_time/3):
    #     #         print "here"
    #     #         break
    #
    #     u.trajectory[0]
    #
    #     print "length of traj = " + str(len(u.trajectory))
    #
    #     num_of_elements = int(round(len(u.trajectory))/3)   # how many elements to have in the list, 1/3 * len(traj)
    #
    #     #num_of_elements = time_counter
    #
    #     list_of_windows = np.zeros(shape=(1, num_of_elements))
    #
    #
    #     print "num of windows used = " + str(num_of_elements)
    #     print list_of_windows[0][1]
    #
    #     for i in range(1, int(num_of_elements)):
    #
    #         try:
    #             print increase
    #         except NameError:
    #             continue
    #
    #         if i <= (num_of_elements / 4):
    #             print "1/4"
    #             increase = 0
    #
    #             if (i + increase) not in list_of_windows:
    #
    #                 element = (i + increase)
    #
    #                 list_of_windows[0][i] = element
    #
    #         elif (num_of_elements / 4) < i <= (num_of_elements / 2):
    #             print "1/4 - 1/2"
    #             increase = 1
    #
    #             if (i + increase) not in list_of_windows:
    #
    #                 element = (i + increase)
    #
    #                 list_of_windows[0][i] = element
    #
    #         elif i > (num_of_elements / 2):
    #             print "1/2 - 1 "
    #             increase = 2
    #
    #             if (i + increase) not in list_of_windows:
    #
    #                 element = (i + increase)
    #
    #                 list_of_windows[0][i] = element
    #
    #     list_of_windows = np.ndarray.tolist(list_of_windows[0])
    #
    #     list_of_windows = filter(lambda a: a != 0.0, list_of_windows)
    #
    #     #print list
    #
    #     window_size_del_t = float(10000) / len(u.trajectory) # ns per frame
    #
    #     print "ns per frame = " + str(window_size_del_t)
    #
    #     return list_of_windows


    # window_list = range(1, 480)

    # window_list = range(1, 700)

    print("input coordinate file: " + str(grofile))
    print("input trajectory file: " + str(trajfile))

    u = mda.Universe(grofile, trajfile)

    max_time = u.trajectory[-1].time

    window_list = increments(u=u)

    protein_sel = u.select_atoms("protein")

    dict_MSD_values = {'MSD_value_dict': {}, 'MSD_std_dict': {}, 'frame_skip_value_list': []}

    dict_individual_MSD_values = {'MSD_value_dict': {}}

    for i, window in enumerate(window_list):

        if window % 10 == 0:
            print("Window = " + str(window))
            stdout.flush()

        counter = 0

        traj_stride_dict = {}

        # individualtrajectory_striding_dictionary = {}



        for ts in u.trajectory[::int(window)]:  # seq[start:end:skip]...

            if counter == 0:  # first frame

                prev_frame_centroid_array_dict = centroid_array_production_protein(protein_sel=protein_sel,
                                                                                   num_protein_copies=options.number_of_proteins)

            else:  # all other frames

                current_frame_centroid_array_dict = centroid_array_production_protein(protein_sel=protein_sel,
                                                                                      num_protein_copies=options.number_of_proteins)

                for particle_name in current_frame_centroid_array_dict.keys():

                    if not particle_name in traj_stride_dict.keys():
                        traj_stride_dict[particle_name] = {'MSD_value_list_centroids': []}

                        # for individual values make dict with tuple of nframes and then zero array
                        # of n_particles -add value of squared distance each iteration and then divide by n frames at end

                        # individualtrajectory_striding_dictionary[particle_name] = {
                        #     'individual_MSD_value_list_centroids': [0, np.zeros(options.number_of_proteins)]}

                    current_delta_array_centroids = prev_frame_centroid_array_dict[particle_name] - \
                                                    current_frame_centroid_array_dict[particle_name]

                    square_darray_centroids = np.square(current_delta_array_centroids)

                    sum_sq_darray_centroids = np.sum(square_darray_centroids, axis=1)

                    # for individual values add here
                    # individualtrajectory_striding_dictionary[particle_name]['individual_MSD_value_list_centroids'][
                    #     0] += 1
                    # individualtrajectory_striding_dictionary[particle_name]['individual_MSD_value_list_centroids'][
                    #     1] += sum_sq_darray_centroids



                    traj_stride_dict[particle_name]['MSD_value_list_centroids'].append(
                        np.average(sum_sq_darray_centroids))  # ensemble avg

                # reset the value of the 'previous' array as you go along:

                prev_frame_centroid_array_dict = current_frame_centroid_array_dict

            counter += 1

        for particle_name, msd_data_subdict in traj_stride_dict.items():  # time avg

            if not particle_name in dict_MSD_values['MSD_value_dict'].keys():  # initialize subdictionaries as needed

                dict_MSD_values['MSD_value_dict'][particle_name] = []

                dict_MSD_values['MSD_std_dict'][particle_name] = []

            dict_MSD_values['MSD_value_dict'][particle_name].append(
                np.average(np.array(traj_stride_dict[particle_name]['MSD_value_list_centroids'])))

            dict_MSD_values['MSD_std_dict'][particle_name].append(
                np.std(np.array(traj_stride_dict[particle_name]['MSD_value_list_centroids'])))

    dict_MSD_values['frame_skip_value_list'].append(window_list)

    # return dict_MSD_values

    # for particle_name, MSD_data_subdictionary in individualtrajectory_striding_dictionary.iteritems():
    #     if not particle_name in dict_individual_MSD_values[
    #         'MSD_value_dict'].keys():  # initialize subdictionaries as needed
    #
    #         dict_individual_MSD_values['MSD_value_dict'][particle_name] = np.zeros(
    #             (len(window_list), options.number_of_proteins))
    #         dict_individual_MSD_values['MSD_value_dict'][particle_name][i] = (
    #         individualtrajectory_striding_dictionary[particle_name]['individual_MSD_value_list_centroids'][1] / float(
    #             individualtrajectory_striding_dictionary[particle_name]['individual_MSD_value_list_centroids'][0]))


    # print dict_individual_MSD_values['MSD_value_dict'][particle_name]

    ### Plotting ###

    fig = plt.figure(figsize=(18, 18))

    ax = fig.add_subplot(111)

    ax.errorbar(window_list, dict_MSD_values['MSD_value_dict'][particle_name],
                yerr=dict_MSD_values['MSD_std_dict'][particle_name], fmt='o', mfc='black', mec='black', markersize=2.0,
                ecolor='skyblue')

    ax.set_title("MSD (" + str(max_time / 1000000) + " $\mu$s)", fontsize=25)
    ax.set_xlabel('Window Size (Frames)', fontsize=20)
    ax.set_ylabel('Mean Squared Displacement ($\AA^2$)', fontsize=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)

    plt.savefig('msd_plot.svg', format='svg', dpi=1000)

    ###############

    return max_time, dict_MSD_values, dict_MSD_values['MSD_value_dict'][particle_name], window_list


def fit_linear(simulation_max_time, time_array, msd_data_array, degrees_of_freedom=2):
    f = open('linear_fit_data.txt', 'w')

    coefficient_dictionary = {1: 2., 2: 4., 3: 6.}  # for mapping degrees_of_freedom to coefficient in fitting equation

    coefficient = coefficient_dictionary[degrees_of_freedom]

    x_data = time_array

    y_data = msd_data_array

    z = np.polyfit(x_data, y_data, 1)

    slope, intercept = z

    # print z
    # print "slope = " + str(slope)
    # print "c = " + str(intercept)

    diff_const = slope / coefficient

    first_half_x_data, second_half_x_data = np.array_split(x_data, 2)
    first_half_y_data, second_half_y_data = np.array_split(y_data, 2)
    slope_first_half, intercept_first_half = np.polyfit(first_half_x_data, first_half_y_data, 1)
    slope_second_half, intercept_second_half = np.polyfit(second_half_x_data, second_half_y_data, 1)

    diffusion_constant_error_estimate = abs(slope_first_half - slope_second_half) / coefficient

    f.write("Diffusion Coefficient = " + str(diff_const) + ", Error = " + str(diffusion_constant_error_estimate) + '\n')
    f.close()

    print "Linear Diffusion Coefficient = " + str(diff_const) + ", Error = " + str(diffusion_constant_error_estimate)

    p = np.poly1d(z)

    fit_x_data = np.linspace(time_array[0], time_array[-1], 100)

    fit_y_data = p(fit_x_data)

    ### Plotting ###

    fig = plt.figure(figsize=(15, 15))

    ax = fig.add_subplot(111)

    if options.dt:
        fit_x_data = (fit_x_data * options.dt) / 1000
        time_array = (time_array * options.dt) / 1000

    ax.plot(fit_x_data, fit_y_data, '-', color='orchid', linewidth=4)
    ax.plot(time_array, msd_data_array, '.', color='black', markersize=12)

    ax.set_title("Linear MSD", fontsize=25)

    if options.dt:
        ax.set_xlabel('$\Delta$t x 10$^3$ (ns)', fontsize=20)
    else:
        ax.set_xlabel('Window Size (Frames)', fontsize=20)

    ax.set_xlabel('Window Size (Frames)', fontsize=20)
    ax.set_ylabel('Mean Squared Displacement / $\AA^2$', fontsize=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)

    plt.savefig('linear_fit_plot.svg', format='svg', dpi=1000)
    # plt.show()

    #################

    return diff_const


def fit_anomalous(simulation_max_time, time_array, msd_data_array, degrees_of_freedom=2):
    g = open('anomalous_fit_data.txt', 'w')

    def fitting_func(time_array, frac_diff_coeff, exponent):
        coeff_dict = {1: 2, 2: 4, 3: 6}
        coefficient = coeff_dict[degrees_of_freedom]

        return coefficient * frac_diff_coeff * (time_array ** exponent)

    optimised_values, covar_params = curve_fit(fitting_func, time_array, msd_data_array)

    # Calc parameters

    stdev_array = np.sqrt(np.diagonal(covar_params))
    fractional_diff_coeff = optimised_values[0]
    stdev_fractional_diff_coeff = stdev_array[0]
    alpha = optimised_values[1]
    stdev_alpha = stdev_array[1]

    g.write("Alpha = " + str(alpha) + ", stdev = " + str(stdev_alpha) + '\n'
                                                                        "Diffusion Coefficient = " + str(
        fractional_diff_coeff) + ", stdev = " + str(stdev_fractional_diff_coeff))
    g.close()

    print "Anomalous Diffusion Coefficient = " + str(fractional_diff_coeff) + \
          ", stdev = " + str(stdev_fractional_diff_coeff)

    # print "alpha from anomalous fit: " + str(alpha)
    # print "alpha stdev: " + str(stdev_alpha)
    # print "frac D from anomalous fit: " + str(fractional_diff_coeff)
    # print "frac D stdev from anomalous fit: " + str (stdev_fractional_diff_coeff)

    fit_x_data = np.linspace(time_array[0], time_array[-1], 100)
    fit_y_data = fitting_func(fit_x_data, *optimised_values)

    ### Plotting ###

    fig = plt.figure(figsize=(15, 15))

    ax = fig.add_subplot(111)

    if options.dt:  # convert to microseconds instead, otherwise keep as "Window Size" (i.e. if ns/frame not known).

        fit_x_data = (fit_x_data * options.dt) / 1000
        time_array = (time_array * options.dt) / 1000

        fit_y_data = fit_y_data / 10  # convert MSD to nm^2
        msd_data_array[:] = [x / 10 for x in msd_data_array]

    ax.plot(fit_x_data, fit_y_data, '-', color='lightseagreen', linewidth=4)
    ax.plot(time_array, msd_data_array, '.', color='black', markersize=12)

    ax.set_title("Anomalous MSD", fontsize=35)

    if options.dt:
        ax.set_xlabel('$\Delta$t ($\mu$s)', fontsize=32)
        ax.set_ylabel('Mean Squared Displacement (nm$^2$)', fontsize=32)
    else:
        ax.set_xlabel('Window Size (Frames)', fontsize=32)
        ax.set_ylabel('Mean Squared Displacement ($\AA^2$)', fontsize=32)

    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=30)

    plt.savefig('anomalous_fit_plot.svg', format='svg', dpi=300)
    # plt.show()

    ################


def fit_log(time_array, msd_data_array, degrees_of_freedom=2):  # To do..

    coeff_dict = {1: 2, 2: 4, 3: 6}

    ln_msd = np.log(coeff_dict[degrees_of_freedom]) + np.log(msd_data_array)

    ln_delt = np.log(time_array)

    # fit_y = (alpha * np.log(time_array)) + np.log(4*D)

    # print fit_y
    #
    # print "diff from log = " + str(D / np.log(coeff_dict[degrees_of_freedom]))
    # print "alpha from log = " + str(alpha)

    fig = plt.figure(figsize=(15, 15))

    ax = fig.add_subplot(111)

    p = ax.plot(np.log(msd_data_array), ln_msd, '.', color='#4d004b')
    # p2 = ax.plot(np.log(time_array), np.log(msd_data_array), '.', color='red')

    # p = ax.plot(ln_msd, ln_delt, '.', color='#4d004b')

    ax.set_title("Log MSD", fontsize=25)

    if options.dt:
        ax.set_xlabel('ln($\Delta$t x 10$^3$) (ns)', fontsize=20)
    else:
        ax.set_xlabel('ln(Window Size (Frames))', fontsize=20)

    ax.set_ylabel('ln(MSD) ($\AA^2$)', fontsize=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    # ax.legend(loc=2)

    plt.savefig('log_plot.svg', format='svg', dpi=1000)


if __name__ == "__main__":
    start_time = time.time()

    parser = argparse.ArgumentParser(description="Calculates MSD at different window sizes & calculates D coefficients")
    parser.add_argument("-c", dest="grofile", type=str, help='the gro file input')
    parser.add_argument("-x", dest="trajfile", type=str,
                        help='the corrected xtc file, processed with nojump option and centered')
    parser.add_argument("-n", dest="number_of_proteins", type=int,
                        help='the number of copies of the protein in the system')
    parser.add_argument("-dt", dest="dt", type=int,
                        help="number of ns for each frame, usually obtianed from gmx trjconv -dt option")
    options = parser.parse_args()

    # sim_length, data_dict, msd_data, del_t_data = msd(grofile='confout_nojump_corrected.gro', trajfile='traj_nojump_corrected.xtc')

    sim_length, data_dict, msd_data, del_t_data = msd(grofile=options.grofile, trajfile=options.trajfile)

    D = fit_linear(simulation_max_time=sim_length, time_array=del_t_data, msd_data_array=msd_data)

    fit_anomalous(simulation_max_time=sim_length, time_array=del_t_data, msd_data_array=msd_data)

    fit_log(time_array=del_t_data, msd_data_array=msd_data)

    print("--- %s seconds ---" % (time.time() - start_time))
