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

    ''' Creates a dictionary containing the "centre" of each protein in the system (i.e. a centroid)'''

    dictionary_centroid_arrays = {}

    full_protein_coord_array = protein_sel.positions

    list_individual_protein_coord_arrays = np.split(full_protein_coord_array, num_protein_copies)

    list_per_protein_centroids = [np.average(protein_coord_array, axis=0) for protein_coord_array in list_individual_protein_coord_arrays]

    dictionary_centroid_arrays['protein'] = np.array(list_per_protein_centroids)

    return dictionary_centroid_arrays # Now have a dictionary of the 'centre' of each protein in the membrane

def msd(grofile, trajfile):


    def make_number_list(Universe):

        num_of_elements = int(round(len(u.trajectory))/8)   # how many elements to have in the list

        list = np.zeros(shape=(1, num_of_elements))

        print "num of elements " + str(num_of_elements)
        print list[0][1]

        for i in range(1, int(num_of_elements)):

            if i not in list:

                if i <= (num_of_elements / 4):

                    increase = 0

                    if (i + increase) not in list:

                        element = (i + increase)

                        list[0][i] = element

                elif (num_of_elements / 4) < i <= ((num_of_elements) / 2):

                    increase = 1

                    if (i + increase) not in list:

                        element = (i + increase)

                        list[0][i] = element

                elif i > ((num_of_elements) / 2):

                    increase = 2

                    if (i + increase) not in list:

                        element = (i + increase)

                        list[0][i] = element

        list = np.ndarray.tolist(list[0])

        list = filter(lambda a: a != 0.0, list)

        print list
        return list


    #window_list = range(1, 480)

    #window_list = range(1, 700)

    u = mda.Universe(grofile, trajfile)

    print u
    print len(u.trajectory)
    max_time = u.trajectory[-1].time

    window_list = make_number_list(Universe=u)

    protein_sel = u.select_atoms("protein")

    dict_MSD_values = {'MSD_value_dict': {}, 'MSD_std_dict': {}, 'frame_skip_value_list': []}


    for window in window_list:

        if window % 10 == 0:

            print "Progress = " + str(window)
            stdout.flush()

        counter = 0

        traj_stride_dict = {}



        for ts in u.trajectory[::int(window)]:  # seq[start:end:skip]...

            if counter == 0:

                prev_frame_centroid_array_dict = centroid_array_production_protein(protein_sel=protein_sel, num_protein_copies=options.number_of_proteins)

            else:

                current_frame_centroid_array_dict = centroid_array_production_protein(protein_sel=protein_sel, num_protein_copies=options.number_of_proteins)

                for particle_name in current_frame_centroid_array_dict.keys():

                    if not particle_name in traj_stride_dict.keys():

                        traj_stride_dict[particle_name] = {'MSD_value_list_centroids': []}

                    current_delta_array_centroids = prev_frame_centroid_array_dict[particle_name] - current_frame_centroid_array_dict[particle_name]

                    square_darray_centroids = np.square(current_delta_array_centroids)

                    sum_sq_darray_centroids = np.sum(square_darray_centroids, axis=1)

                    traj_stride_dict[particle_name]['MSD_value_list_centroids'].append(np.average(sum_sq_darray_centroids)) #ensemble avg

                prev_frame_centroid_array_dict = current_frame_centroid_array_dict

            counter += 1

        for particle_name, msd_data_subdict in traj_stride_dict.items(): # time avg

            if not particle_name in dict_MSD_values['MSD_value_dict'].keys():  # initialize subdictionaries as needed

                dict_MSD_values['MSD_value_dict'][particle_name] = []

                dict_MSD_values['MSD_std_dict'][particle_name] = []

            dict_MSD_values['MSD_value_dict'][particle_name].append(np.average(np.array(traj_stride_dict[particle_name]['MSD_value_list_centroids'])))

            dict_MSD_values['MSD_std_dict'][particle_name].append(np.std(np.array(traj_stride_dict[particle_name]['MSD_value_list_centroids'])))

    dict_MSD_values['frame_skip_value_list'].append(window_list)

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

    ax.plot(fit_x_data, fit_y_data, '-', color='orchid', linewidth=3)
    ax.plot(time_array, msd_data_array, '.', color='black')

    ax.set_title("Linear MSD Fit (" + str(simulation_max_time /1000000) + " $\mu$s)", fontsize=25)
    ax.set_xlabel('Window Size (Frames)', fontsize=20)
    ax.set_ylabel('Mean Squared Displacement / $\AA^2$', fontsize=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)

    plt.savefig('linear_fit_plot.svg', format='svg', dpi=1000)
    #plt.show()

    #################

    return diff_const

def fit_anomalous(simulation_max_time, time_array, msd_data_array, degrees_of_freedom=2):

    g = open('anomalous_fit_data.txt', 'w')

    def fitting_func(time_array, frac_diff_coeff, exponent):
        coeff_dict = {1:2, 2:4, 3:6}
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
            "Diffusion Coefficient = " + str(fractional_diff_coeff) + ", stdev = " + str(stdev_fractional_diff_coeff))
    g.close()

    print "Anomalous Diffusion Coefficient = " + str(fractional_diff_coeff) +\
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

    ax.plot(fit_x_data, fit_y_data, '-', color='lightseagreen', linewidth=3)
    ax.plot(time_array, msd_data_array, '.', color='black')

    ax.set_title("Anomalous MSD Fit (" + str(simulation_max_time / 1000000) + " $\mu$s)", fontsize=25)
    ax.set_xlabel('Window Size (Frames)', fontsize=20)
    ax.set_ylabel('Mean Squared Displacement / $\AA^2$', fontsize=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)

    plt.savefig('anomalous_fit_plot.svg', format='svg', dpi=1000)
    #plt.show()

    ################

def fit_log(time_array, msd_data_array, degrees_of_freedom=2): # To do..

    coeff_dict = {1: 2, 2: 4, 3: 6}

    ln_msd = np.log(coeff_dict[degrees_of_freedom]) + np.log(msd_data_array)

    fig = plt.figure()

    ax = fig.add_subplot(111)

    p = ax.plot(np.log(msd_data_array), ln_msd, '-', color='black')
    #p2 = ax.plot(time_array, msd_data_array, '.')

    c = ax.set_xlabel('ln(Window Size (Frames))')
    c = ax.set_ylabel('ln(MSD) / $\AA^2$')
    c = ax.legend(loc=2)

    plt.show()

if __name__ == "__main__":

    start_time = time.time()

    parser = argparse.ArgumentParser(description="Calculates MSD at different window sizes & calculates D coefficients")
    parser.add_argument("-c", dest="grofile", type=str, help='the gro file input')
    parser.add_argument("-x", dest="trajfile", type=str, help='the corrected xtc file, processed with nojump option and centered')
    parser.add_argument("-n", dest="number_of_proteins", type=int, help='the number of copies of the protein in the system')
    options = parser.parse_args()

    #sim_length, data_dict, msd_data, del_t_data = msd(grofile='confout_nojump_corrected.gro', trajfile='traj_nojump_corrected.xtc')
    sim_length, data_dict, msd_data, del_t_data = msd(grofile=options.grofile, trajfile=options.trajfile)

    D = fit_linear(simulation_max_time=sim_length, time_array=del_t_data, msd_data_array=msd_data)

    fit_anomalous(simulation_max_time=sim_length, time_array=del_t_data, msd_data_array=msd_data)

    #fit_log(time_array=del_t_data, msd_data_array=msd_data)

    print("--- %s seconds ---" % (time.time() - start_time))
