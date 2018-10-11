import MDAnalysis as mda
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
import argparse

scale_factor = 20

def get_universe(gro_file, traj_file):

    u = mda.Universe(gro_file, traj_file)

    return u


def get_principal_axes(universe, selection):

    # Methodology taken from Beckstein
    # https://stackoverflow.com/questions/49239475/
    # how-to-use-mdanalysis-to-principal-axes-and-moment-of-inertia-with-a-group-of-at/49268247#49268247

    CA = universe.select_atoms(selection)

    I = CA.moment_of_inertia()
    # UT = CA.principal_axes()
    # U = UT.T

    values, evecs = np.linalg.eigh(I)
    indices = np.argsort(values)
    U = evecs[:, indices]

    # Lambda = U.T.dot(I.dot(U))
    #
    # print(Lambda)
    # print(np.allclose(Lambda - np.diag(np.diagonal(Lambda)), 0))

    # TODO: Need to ensure that the z vector of the system is orientated correctly. (not exactly a prob with code)

    return U


def get_com(universe, selection):

    com = universe.select_atoms(selection).center_of_mass()

    return com

#
# def correction(v):
#     '''
#     :param v: this is a principal axis vector
#     :return: a "corrected vector", this will vary depending on what vector you are looking at (pitch, roll or yaw)
#             for example: if looking at pitch, this will return a vector that is +ve or -ve in z and always +ve in
#             x & y.
#     '''
#
#     v = np.array(v)
#
#     if v[0] < 0 and v[1] < 0: # if x & y -ve
#         #print("\nvector before correction: ", v)
#
#         #v[2] = v[2] * -1
#         v = v * -1
#         #print("vector after correction: ", v)
#
#
#     if v[0] < 0 and v[1] > 0: # if x -ve but y +ve
#         #v[2] = v[2] * -1
#         v = v * -1
#
#     return v

class correction(object):

    '''
    :param v: this is a principal axis vector
    :return: a "corrected vector", this will vary depending on what vector you are looking at (pitch, roll or yaw)
            for example: if looking at pitch, this will return a vector that is +ve or -ve in z and always +ve in
            x & y.
    '''

    def __init__(self, v):

        self.v = v
        #print(self.v)


    def pitch(self):

        if self.v[0] < 0 and self.v[1] < 0: # if x & y -ve


            v_corr = self.v * -1


        if self.v[0] < 0 and self.v[1] > 0: # if x -ve but y +ve
            #v[2] = v[2] * -1
            v_corr = self.v * -1

        else:
            v_corr = self.v


        return v_corr


    def roll(self):

        if self.v[1] < 0 and self.v[2] < 0:  # if y & z -ve

            v_corr = self.v * -1


        if self.v[2] < 0 and self.v[1] > 0:  # if z -ve but y +ve

            v_corr = self.v * -1

        else:
            v_corr = self.v

        return v_corr


def dir_angle(v1, v2, prev_angle, i, measure):

    # Since the vectors between frames can change direction (due to the calc of the PA's) then we add a cut off
    # the cut off depends on what measurement we are taking. larger cut offs for pitch (as they are slower)
    # and smaller cut offs for roll since they can be faster (i.e. a protein can roll faster than it can change a large
    # amount in pitch) - this is a bit of a hack.
    # TODO: find a way to counteract the change in vector direction w/o arbitrary cut offs.

    if measure == "pitch":
        cut_off = 30
    elif measure == "roll":
        cut_off = 30
    else:
        cut_off = 30

    print(cut_off)

    if i == 0:

        direction_cosine = np.arccos(np.dot(v1, v2)) / np.linalg.norm(v1)
        return direction_cosine



    # if np.dot(v1, v2) < 0:
    #     print ("obtuse!", 180*angle/np.pi)

    # if acute:
    #     return angle
    # else:
    #     return (2 * np.pi) - angle

    # See: https://stackoverflow.com/questions/39497496/angle-between-two-vectors-3d-python

    # See: https: // en.wikipedia.org / wiki / Direction_cosine
    #direction_cosine = np.arccos(np.abs(np.dot(v1, v2))) / np.linalg.norm(v1) # * np.linalg.norm(v2)


    elif i > 0:

        #angle_prev = np.arccos(np.dot(axis_prev, v2)) / np.linalg.norm(axis_prev)

        direction_cosine = np.arccos(np.dot(v1, v2)) / np.linalg.norm(v1)  # * np.linalg.norm(v2)

        # to correct for if the PA calc orientates the vector the other way around, need to check if this is a genuine rot
        if measure == "pitch":
            print("Current angle = ", np.rad2deg(direction_cosine))
            print("Previous angle = ", np.rad2deg(prev_angle))
            print("Difference between prev angle: ", np.rad2deg(direction_cosine) - np.rad2deg(prev_angle))

        #if np.abs(direction_cosine - prev_angle) >= np.deg2rad(90):

        #if np.abs(direction_cosine - prev_angle) >= np.deg2rad(cut_off):
        if direction_cosine - prev_angle >= np.deg2rad(cut_off):

            if measure == "pitch":
                print("")
                print("orig angle:", np.rad2deg(direction_cosine))

            direction_cosine_corrected = np.deg2rad(180) - direction_cosine

            if measure == "pitch":
                print("new angle:", np.rad2deg(direction_cosine_corrected))
                print("")

            return direction_cosine_corrected

        else:

            return direction_cosine

def make_direction_cosine_matrix(axes_set):

    ex = [1, 0, 0]
    ey = [0, 1, 0]
    ez = [0, 0, 1]

    u = axes_set[:, 0] # 1st PA
    v = axes_set[:, 1] # 2nd PA
    w = axes_set[:, 2] # 3rd PA

    u_alpha = dir_angle(u, ex)
    v_alpha = dir_angle(v, ex)
    w_alpha = dir_angle(w, ex)

    u_beta = dir_angle(u, ey)
    v_beta = dir_angle(v, ey)
    w_beta = dir_angle(w, ey)

    u_gamma = dir_angle(u, ez)
    v_gamma = dir_angle(v, ez)
    w_gamma = dir_angle(w, ez)

    c_matrix = np.matrix([[u_alpha, v_alpha, w_alpha],
                          [u_beta, v_beta, w_beta],
                          [u_gamma, v_gamma, w_gamma]])

    #print(c_matrix)

    return c_matrix


def vis_axes(vis, axes_data, center):

    # PAs returned in row format. i.e. the PAs are the columns of axes. To get the 1st PA: axes[:, 0] etc
    axis1 = axes_data[:, 0]
    axis2 = axes_data[:, 1]
    axis3 = axes_data[:, 2]

    # print("")
    # print("VMD")
    # print(axis1)
    #axis1 = correction(axis1)
    axis1 = correction(axis1).pitch()
    axis2 = correction(axis2).roll()

    print(" AXIS 1:", axis1)
    print("AXIS 2:", axis2 )
    # this is a visual fix, the z component may actually point down as this is from the PA calc.
    # this just ensures that the the vector doesn't fly around in the trajectory when visualising it.
    # if axis1[2] < 0:
    #     axis1[2] = axis1[2] * -1

    if axis1[2] < 0:
        axis1 = axis1 * -1

    if axis2[0] < 0:
        axis2 = axis2 * -1


    # this is a visual fix, the x component may actually point down as this is from the PA calc.
    # this just ensures that the the vector doesn't fly around in the trajectory when visualising it.
    # if axis2[0] < 0:
    #     axis2[0] = axis2[1] * -1



    if vis == 'vmd':

        i, j, k = 1, 1, 1

        output = open('pa_vectors.pdb', 'a')

        for i in range(0, (3 * scale_factor)):
            tmp = "ATOM    {0:3d}  CA  ALA A {1:3d}    {2:8.3f}{3:8.3f}{4:8.3f}  1.00  0.00\n".format(i, i, center[0] + (axis1[0] * i),
                                                                                                    center[1] + (axis1[1] * i),
                                                                                                    center[2] + (axis1[2] * i))
            output.write(tmp)
            #print(tmp, file=output)
            i += 1
        output.write("TER\n")
        #print("TER", file=output)

        for j in range(0, (2 * scale_factor)):
            tmp2 = "ATOM    {0:3d}  CA  ALA B {1:3d}    {2:8.3f}{3:8.3f}{4:8.3f}  1.00  0.00\n".format(j, j, center[0] + (axis2[0] * j),
                                                                                                    center[1] + (axis2[1] * j),
                                                                                                    center[2] + (axis2[2] * j))
            output.write(tmp2)
            #print(tmp2, file=output)
            j += 1
        output.write("TER\n")
        #print("TER", file=output)

        for k in range(0, (1 * scale_factor)):
            tmp3 = "ATOM    {0:3d}  CA  ALA C {1:3d}    {2:8.3f}{3:8.3f}{4:8.3f}  1.00  0.00\n".format(k, k, center[0] + (axis3[0] * k),
                                                                                                    center[1] + (axis3[1] * k),
                                                                                                    center[2] + (axis3[2] * k))
            output.write(tmp3)
            #print(tmp3, file=output)
           # k += 1
        output.write("TER\nENDMDL\n")
        #print("TER\nENDMDL", file=output)

        output.close()

    elif vis == 'pymol':

        # --------------------------------------------------------------------------
        # center axes to the geometric center of the molecule
        # and rescale them by order of eigen values
        # --------------------------------------------------------------------------

        # the large vector is the first principal axis
        point1 = 3 * scale_factor * axis1 + center
        # the medium vector is the second principal axis
        point2 = 2 * scale_factor * axis2 + center
        # the small vector is the third principal axis
        point3 = 1 * scale_factor * axis3 + center

        #pymol_name = pdb_name.replace(".pdb", "_axes.pml")

        pymol_name = ("test_axes.pml")
        with open(pymol_name, "w") as pymol_file:
            pymol_file.write(
                """
                from cgo import *
                axis1=  [ BEGIN, LINES, COLOR, 1.0, 0.0, 0.0, \
                VERTEX, %8.3f, %8.3f, %8.3f, VERTEX, %8.3f, %8.3f, %8.3f, END ]
                axis2=  [ BEGIN, LINES, COLOR, 0.0, 1.0, 0.0, \
                VERTEX, %8.3f, %8.3f, %8.3f, VERTEX, %8.3f, %8.3f, %8.3f, END ]
                axis3=  [ BEGIN, LINES, COLOR, 0.0, 0.0, 1.0, \
                VERTEX, %8.3f, %8.3f, %8.3f, VERTEX, %8.3f, %8.3f, %8.3f, END ]
                cmd.load_cgo(axis1, 'axis1')
                cmd.load_cgo(axis2, 'axis2')
                cmd.load_cgo(axis3, 'axis3')
                cmd.set('cgo_line_width', 4)
                """ % ( \
                    center[0], center[1], center[2], point1[0], point1[1], point1[2], \
                    center[0], center[1], center[2], point2[0], point2[1], point2[2], \
                    center[0], center[1], center[2], point3[0], point3[1], point3[2]))

        # --------------------------------------------------------------------------
        # create .pml script for nice rendering in Pymol
        # output usage
        # --------------------------------------------------------------------------
        print("\nFirst principal axis (in red)")
        # print("coordinates: ", axis1)
        # print("eigen value: ", eval1)

        print("\nSecond principal axis (in green)")
        # print("coordinates:", axis2)
        # print("eigen value:", eval2)

        print("\nThird principal axis (in blue)")

def run_traj(u, chain_info, skip_val):

    k = 0

    sel_resids = ' or resid '.join(map(str, chain_info['chain 0']['resids']))

    sel = "name CA and (resid " + str(sel_resids) + ")"

    #print("COM", get_com(u, sel))

    # create coordinates array
    coord = np.array(u.select_atoms(sel).atoms.positions, float)

    # # compute geometric center
    # center = np.mean(coord, 0)
    # #print("Coordinates of the geometric center:\n", center)
    #
    # # center with geometric center
    # coord = coord - get_com(u, sel)

    start_axes = get_principal_axes(u, sel)

    prev_axis1 = start_axes[:, 0]
    prev_axis2 = start_axes[:, 1]
    # axis3 = axes[:, 2]

    prev_axis1 = correction(prev_axis1).pitch()
    prev_axis2 = correction(prev_axis2).roll()

    prev_angle1 = dir_angle(prev_axis1, [0, 0, 1], prev_angle=None, i=k, measure="pitch")
    prev_angle2 = dir_angle(prev_axis2, [1, 0, 0], prev_angle=None, i=k, measure="roll")

    for ts in u.trajectory[::skip_val]:

        print("Frame = ", ts.frame, ", time = ", ts.time / 1000)
        print("")

        # At each frame calculate the angle that each chain has moved through
        for chain in chain_info:

            sel_resids = ' or resid '.join(map(str, chain_info[chain]['resids']))

            sel = "name CA and (resid " + str(sel_resids) + ")"

            pa_array = get_principal_axes(u, sel)

            #dir_cosine_matrix = make_direction_cosine_matrix(pa_array)


            # PAs returned in row format. i.e. the PAs are the columns of axes. To get the 1st PA: axes[:, 0] etc

            ax1 = pa_array[:, 0]  # The principal axis
            #
            # #ax1 = correction(ax1)
            ax1 = correction(ax1).pitch()
            #
            # print("")
            # print("ax1 = ", ax1)
            # print("")

            ax2 = pa_array[:, 1]  # The 2nd principal axis

            ax2 = correction(ax2).roll()
            #
            # print("")
            # print("ax2 = ", ax2)
            # print("")

            #ax2 = pa_array[:, 1]  # The 2nd principal axis
            #ax3 = pa_array[:, 2]  # The 3rd principal axis


            if options.vector_traj == True:

                # create coordinates array
                coord = np.array(u.select_atoms(sel).atoms.positions, float)

                # compute geometric center
                center = np.mean(coord, 0)

                vis_axes(vis='vmd', axes_data=pa_array, center=center)


            # measuring to standard z, y, x vectors. We define x as the principal axis, y as 2nd and z as third (roll, pitch, yaw).

            pa_1_measured_angle = dir_angle(ax1, [0, 0, 1], prev_angle=prev_angle1, i=k, measure="pitch") # get the angle b/w PA1 and the z axis (pitch)
            pa_2_measured_angle = dir_angle(ax2, [1, 0, 0], prev_angle=prev_angle2, i=k, measure="roll") # get the angle b/w PA2 and the x axis (roll)


            #pa_2_measured_angle = dir_angle(ax2, [0, 1, 0]) # get the angle between PA2 and the y axis (yaw)
            #pa_3_measured_angle = dir_angle(ax3, [1, 0, 0]) # get the angle between PA3 and the x axis (roll)

            # pa_1_measured_angle_2 = dir_cosine_matrix[2, 0]  # get the angle between PA1 and the z axis (pitch)
            # pa_2_measured_angle_2 = dir_cosine_matrix[1, 1]   # get the angle between PA2 and the y axis (yaw)
            # pa_3_measured_angle_2 = dir_cosine_matrix[0, 2]  # get the angle between PA3 and the x axis (roll)

           # print("dir cosine PA1 = ", pa_1_measured_angle)
            print("dir cosine PA2 = ", np.rad2deg(pa_2_measured_angle))

            chain_info[chain]['angle_pa1'].append((ts.time / 1000, np.rad2deg(pa_1_measured_angle)))
            chain_info[chain]['angle_pa2'].append((ts.time / 1000, np.rad2deg(pa_2_measured_angle)))
            #chain_info[chain]['angle_pa3'].append((ts.time / 1000, np.rad2deg(pa_3_measured_angle)))

            # chain_info[chain]['angle_pa1'].append((ts.time / 1000, np.rad2deg(pa_1_measured_angle_2)))
            # chain_info[chain]['angle_pa2'].append((ts.time / 1000, np.rad2deg(pa_2_measured_angle_2)))
            # chain_info[chain]['angle_pa3'].append((ts.time / 1000, np.rad2deg(pa_3_measured_angle_2)))

            prev_angle1 = pa_1_measured_angle
            prev_angle2 = pa_2_measured_angle

        k += 1
    #print(chain_info)

        # if ts.time / 1000 == 200.0:
        #
        #     plt.plot(a, label='a')
        #     plt.plot(b, label='b')
        #     plt.plot(c, label='c')
        #     plt.legend()
        #     plt.show()

    return pa_array, chain_info


def read_stride(file, num_prots, chain_length, specific_region):
    ''' This reads the output from running stride structure.pdb, this is used to identify beta sheets and alpha
    helices. Due to flexible loops the calculated principal axes can differ so using more stable regions can give
    less noisy data.

    Prior to this function the user must run "stride file.pdb > stride_file.txt"

    '''

    x = []
    ss_list = ['E']#, 'G', 'H', 'I']
    with open(file) as f:

        for line in f.readlines():

            if line.splitlines()[0][0] == 'A':

                if line.splitlines()[0][24] in ss_list:
                    res = (line.splitlines()[0][12], line.splitlines()[0][13], line.splitlines()[0][14])
                    x.append(int(''.join(res)))

    #chain_labels = ['A', 'B', 'C']#, 'D', 'E']

    chain_dict = {'chain ' + str(i): {'resids': [],
                                      'angle_pa1': [],
                                      'angle_pa2': [],
                                      'angle_pa3': [],
                                      } for i in range(num_prots)}

    print (x)
    print(len(x))
    print("HERE")
    print(num_prots, chain_length)

    print(chain_dict)

    for i in range(num_prots):

        if i == 0:

            chain_dict['chain ' + str(i)]['resids'] = [t for t in x if 1 <= t <= chain_length]

        # Need to test below works on a multi chain system
        else:

            chain_dict['chain ' + str(i)]['resids'] = [t for t in x if
                                               (i * chain_length) + 1 <= t <= ((i+1) * chain_length)]  # 124 to 246


    #chain_dict['chain 1']['resids'] = [t for t in x if 1 <= t <= chain_length] # 1 to 123

    #print(chain_dict)

    #chain_dict['chain 2']['resids'] = [t for t in x if (chain_length + 1) <= t <= (2 * chain_length)] # 124 to 246

    #chain_dict['chain 3']['resids'] = [t for t in x if (1 + (2 * chain_length)) <= t <= (3 * chain_length)] # 247 to 369

    print(chain_dict)

    return chain_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculates the orientation of a user defined region of a protein")
    parser.add_argument("-c", dest="gro_file", help='the coordinate file [.gro]')
    parser.add_argument("-f", dest="xtc_file", help='a corrected trajectory file, pbc artifacts removed and the protein centered')

    parser.add_argument("-sel", dest="selection", type=str, help='the range of resids to use, in the form of A:B where A and B are integers')

    parser.add_argument("-n", dest="num_of_proteins", type=int, help='the number of protein copies in the system')

    parser.add_argument("-s", dest="skip", type=int, help="the number of frames to skip", default=1)

    parser.add_argument("-vtraj", dest="vector_traj", type=bool, help="set to True if you want a trajectory of the vectors")

    #parser.add_argument("-pt", dest="plot_type", type=str, help="type of plot, either 'single' or 'multiple'")

    options = parser.parse_args()

start_of_protein = int(options.selection.split()[-1].split(':')[0])
end_of_protein = int(options.selection.split()[-1].split(':')[1])

universe = get_universe(options.gro_file, options.xtc_file)

#prot_length = len(universe.select_atoms("protein and name CA"))
prot_length = len(universe.select_atoms("resid " + str(start_of_protein) + ":" + str(end_of_protein) + " and name CA"))

#chain_dict = read_stride('stride_file.txt', options.num_of_proteins, int(prot_length / options.num_of_proteins), False)

chain_dict = read_stride('stride_file.txt', options.num_of_proteins, int(prot_length), False)

princ_axes, data = run_traj(universe, chain_dict, options.skip)








### plotting (make into function) ###

my_cmap = ListedColormap(sns.color_palette())

if options.num_of_proteins > 1:

    fig, ax = plt.subplots(len(data), sharex=True, sharey=True)


    for i, chain in enumerate(data):

        # ax[i].scatter(*zip(*data[chain]['angle']), label=str(chain), c=col[i])
        # ax[i].legend()

        ax[i].scatter(*zip(*data[chain]['angle_pa1']), label=str(chain), c=my_cmap, alpha=0.8)
        ax[i].scatter(*zip(*data[chain]['angle_pa2']), label=str(chain), c=my_cmap, alpha=0.8)
        ax[i].scatter(*zip(*data[chain]['angle_pa3']), label=str(chain), c=my_cmap, alpha=0.8)
        ax[i].legend()

    #plt.scatter(*zip(*dat))

    # Add titles, legend, and save
    #ax.set(title=" ", xlabel="Time (ns)", ylabel="Tilt Angle (" + r'$^{\circ}$' + ")")

    fig.text(0.5, 0.04, 'Time (ns)', ha='center')
    fig.text(0.04, 0.5, "Angle (" + r'$^{\circ}$' + ")", va='center', rotation='vertical')

    plt.savefig('protein_tilt_pa_to_z_all_chains.svg', format='svg')

    #plt.show()

else:

    fig, ax = plt.subplots()

    # set the colour palette based on the Seabron colorblind colours
    ax.set_prop_cycle('color', sns.color_palette("colorblind", 3))

    ax.plot(*zip(*data['chain 0']['angle_pa1']), alpha=0.8, label='princ axis 1 (pitch)')
    ax.plot(*zip(*data['chain 0']['angle_pa2']), alpha=0.8, label='princ axis 2 (roll)')
    ax.plot(*zip(*data['chain 0']['angle_pa3']), alpha=0.8, label='princ axis 3 (yaw)')
    ax.legend()

    fig.text(0.5, 0.04, 'Time (ns)', ha='center')
    fig.text(0.04, 0.5, "Angle (" + r'$^{\circ}$' + ")", va='center', rotation='vertical')

    plt.savefig('protein_tilt_pa_to_z.svg', format='svg')

    #plt.show()

    # fig, ax = plt.subplots()
    #
    # # set the colour palette based on the Seabron colorblind colours
    # ax.set_prop_cycle('color', sns.color_palette("colorblind", 3))
    #
    # ax.scatter(*zip(*data['chain 0']['angle_pa1']), alpha=0.8, label='princ axis 1 (pitch)')
    # ax.scatter(*zip(*data['chain 0']['angle_pa2']), alpha=0.8, label='princ axis 2 (yaw)')
    # ax.scatter(*zip(*data['chain 0']['angle_pa3']), alpha=0.8, label='princ axis 3 (roll)')
    # ax.legend()
    #
    # fig.text(0.5, 0.04, 'Time (ns)', ha='center')
    # fig.text(0.04, 0.5, "Angle (" + r'$^{\circ}$' + ")", va='center', rotation='vertical')
    #
    # plt.savefig('protein_tilt_pa_to_z_2.svg', format='svg')













# print("COM", get_com(u, sel))
#
# #create coordinates array
# coord = np.array(u.select_atoms(sel).atoms.positions, float)
#
# # compute geometric center
# center = np.mean(coord, 0)
# print("Coordinates of the geometric center:\n", center)
#
# # center with geometric center
# coord = coord - get_com(u, sel)


# def check_argument(arguments):
#     """
#     Check if filename passed as argument exists.
#     Parameters
#     ----------
#     arguments : list
#         list of arguments passed to the script
#     Returns
#     -------
#     string
#         file name
#     """
#     if len(arguments) == 3:
#         file_name = arguments[1]
#
#     elif len(arguments) == 2:
#         file_name = arguments[1]
#
#     else:
#         message = """
#         ERROR: missing pdb filename as argument
#         usage: %s file.pdb""" %(arguments[0])
#         sys.exit(message)
#
#     # check if argument is an existing file
#     if not os.path.exists(file_name):
#         sys.exit("ERROR: file %s does not seem to exist" %(file_name))
#
#     return file_name

# print(sys.argv)
# pdb_name = check_argument(sys.argv)

# axes = princ_axes
# print(axes)

# axes = get_principal_axes(u, sel)
#
# axis1 = axes[:,0] # ?
# axis2 = axes[:,1] # ?
# axis3 = axes[:,2]
#
# # axis1 = ref_princ_axes[0]
# # axis2 = ref_princ_axes[1]
# # axis3 = ref_princ_axes[2]
#
# print(axis1, axis2, axis3)
#
# #--------------------------------------------------------------------------
# # center axes to the geometric center of the molecule
# # and rescale them by order of eigen values
# #--------------------------------------------------------------------------
# # the large vector is the first principal axis
# point1 = 3 * scale_factor * axis1 + center
# # the medium vector is the second principal axis
# point2 = 2 * scale_factor * axis2 + center
# # the small vector is the third principal axis
# point3 = 1 * scale_factor * axis3 + center
#
# #--------------------------------------------------------------------------
# # create .pml script for a nice rendering in Pymol
# #--------------------------------------------------------------------------
# #pymol_name = pdb_name.replace(".pdb", "_axes.pml")
# pymol_name = pdb_name.replace(".pdb", "_axes.pml")
# with open(pymol_name, "w") as pymol_file:
#     pymol_file.write(
#         """
#         from cgo import *
#         axis1=  [ BEGIN, LINES, COLOR, 1.0, 0.0, 0.0, \
#         VERTEX, %8.3f, %8.3f, %8.3f, VERTEX, %8.3f, %8.3f, %8.3f, END ]
#         axis2=  [ BEGIN, LINES, COLOR, 0.0, 1.0, 0.0, \
#         VERTEX, %8.3f, %8.3f, %8.3f, VERTEX, %8.3f, %8.3f, %8.3f, END ]
#         axis3=  [ BEGIN, LINES, COLOR, 0.0, 0.0, 1.0, \
#         VERTEX, %8.3f, %8.3f, %8.3f, VERTEX, %8.3f, %8.3f, %8.3f, END ]
#         cmd.load_cgo(axis1, 'axis1')
#         cmd.load_cgo(axis2, 'axis2')
#         cmd.load_cgo(axis3, 'axis3')
#         cmd.set('cgo_line_width', 4)
#         """ %( \
#                 center[0], center[1], center[2], point1[0], point1[1], point1[2], \
#                 center[0], center[1], center[2], point2[0], point2[1], point2[2], \
#                 center[0], center[1], center[2], point3[0], point3[1], point3[2]))
#
# #--------------------------------------------------------------------------
# # create .pml script for nice rendering in Pymol
# # output usage
# #--------------------------------------------------------------------------
# print("\nFirst principal axis (in red)")
# # print("coordinates: ", axis1)
# # print("eigen value: ", eval1)
#
# print("\nSecond principal axis (in green)")
# # print("coordinates:", axis2)
# # print("eigen value:", eval2)
#
# print("\nThird principal axis (in blue)")
# print("coordinates:", axis3)
# print("eigen value:", eval3)

#     print("\nYou can view principal axes with PyMOL:")
# print("pymol %s %s" %(pymol_name, pdb_name))











# pdb_name = check_argument(sys.argv)
#
# axes_2 = ref_princ_axes
# print(axes)
#
# refaxis1 = axes_2[:,0] # ?
# refaxis2 = axes_2[:,1] # ?
# refaxis3 = axes_2[:,2]
#
# # axis1 = ref_princ_axes[0]
# # axis2 = ref_princ_axes[1]
# # axis3 = ref_princ_axes[2]
#
# print(refaxis1, refaxis2, refaxis3)
#
# #--------------------------------------------------------------------------
# # center axes to the geometric center of the molecule
# # and rescale them by order of eigen values
# #--------------------------------------------------------------------------
# # the large vector is the first principal axis
# point1 = 3 * scale_factor * refaxis1 + center
# # the medium vector is the second principal axis
# point2 = 2 * scale_factor * refaxis2 + center
# # the small vector is the third principal axis
# point3 = 1 * scale_factor * refaxis3 + center
#
# #--------------------------------------------------------------------------
# # create .pml script for a nice rendering in Pymol
# #--------------------------------------------------------------------------
# #pymol_name = pdb_name.replace(".pdb", "_axes.pml")
# pymol_name = pdb_name.replace(".pdb", "_axes_ref.pml")
# with open(pymol_name, "w") as pymol_file:
#     pymol_file.write(
#         """
#         from cgo import *
#         refaxis1=  [ BEGIN, LINES, COLOR, 1.0, 0.0, 0.0, \
#         VERTEX, %8.3f, %8.3f, %8.3f, VERTEX, %8.3f, %8.3f, %8.3f, END ]
#         refaxis2=  [ BEGIN, LINES, COLOR, 0.0, 1.0, 0.0, \
#         VERTEX, %8.3f, %8.3f, %8.3f, VERTEX, %8.3f, %8.3f, %8.3f, END ]
#         refaxis3=  [ BEGIN, LINES, COLOR, 0.0, 0.0, 1.0, \
#         VERTEX, %8.3f, %8.3f, %8.3f, VERTEX, %8.3f, %8.3f, %8.3f, END ]
#         cmd.load_cgo(refaxis1, 'refaxis1')
#         cmd.load_cgo(refaxis2, 'refaxis2')
#         cmd.load_cgo(refaxis3, 'refaxis3')
#         cmd.set('cgo_line_width', 4)
#         """ %( \
#                 center[0], center[1], center[2], point1[0], point1[1], point1[2], \
#                 center[0], center[1], center[2], point2[0], point2[1], point2[2], \
#                 center[0], center[1], center[2], point3[0], point3[1], point3[2]))
#
# #--------------------------------------------------------------------------
# # create .pml script for nice rendering in Pymol
# # output usage
# #--------------------------------------------------------------------------
# print("\nFirst principal axis (in red)")
# # print("coordinates: ", axis1)
# # print("eigen value: ", eval1)
#
# print("\nSecond principal axis (in green)")
# # print("coordinates:", axis2)
# # print("eigen value:", eval2)
#
# print("\nThird principal axis (in blue)")
# # print("coordinates:", axis3)
# # print("eigen value:", eval3)
#
# #     print("\nYou can view principal axes with PyMOL:")
# print("pymol %s %s" %(pymol_name, pdb_name))