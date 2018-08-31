import MDAnalysis as mda
import numpy as np
import sys
import os
import matplotlib.pyplot as plt

scale_factor = 20

def get_universe(gro_file, traj_file):

    u = mda.Universe(gro_file, traj_file)

    return u


def get_principal_axes(universe, selection):

    # https://stackoverflow.com/questions/49239475/
    # how-to-use-mdanalysis-to-principal-axes-and-moment-of-inertia-with-a-group-of-at/49268247#49268247

    CA = universe.select_atoms(selection)

    I = CA.moment_of_inertia()
    #UT = CA.principal_axes()
    #U = UT.T

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


def angle(v1, v2, acute):

    # See: https://stackoverflow.com/questions/39497496/angle-between-two-vectors-3d-python

    angle = np.arccos(np.dot(v1, v2)) / np.linalg.norm(v1) * np.linalg.norm(v2)

    if acute:
        return angle
    else:
        return (2 * np.pi) - angle


def run_traj(u, chain_info):

    for ts in u.trajectory[::10]:

        print(ts.time / 1000)

        # At each frame calcualte the angle that each chain has moved through
        for chain in chain_info:

            sel_resids = ' or resid '.join(map(str, chain_info[chain]['resids']))

            sel = "name CA and (resid " + str(sel_resids) + ")"

            #temp

            if chain == 'Chain A':
                sel = "name CA and resid 1:123"

            elif chain == 'Chain B':
                sel = "name CA and resid 124:246"

            elif chain == 'Chain C':
                sel = "name CA and resid 247:369"

            else:
                break

            pa_array = get_principal_axes(u, sel)

            ax1 = pa_array[0]# The principal axis

            measured_angle = angle(ax1, [0, 0, 1], True)

            if ax1[2] < 0: # To check if the vector is pointing 180 deg, if it is then correct
                measured_angle = np.pi - measured_angle

            #print("Angle = ", (180 * measured_angle) / np.pi)

            chain_info[chain]['angle'].append((ts.time / 1000, 180 * measured_angle / np.pi))

    return pa_array, chain_info


def read_stride(file, num_prots, chain_length, specific_region):
    ''' This reads the output from running stride structure.pdb, this is used to identify beta sheets and alpha
    helices. Due to flexible loops the calculated principal axes can differ so using more sable regions can give
    less noisy data.'''

    x = []
    ss_list = ['E', 'G', 'H', 'I']
    with open(file) as f:

        for line in f.readlines():

            if line.splitlines()[0][0] == 'A':

                if line.splitlines()[0][24] in ss_list:
                    res = (line.splitlines()[0][12], line.splitlines()[0][13], line.splitlines()[0][14])
                    x.append(int(''.join(res)))

    chain_labels = ['A', 'B', 'C']#, 'D', 'E']

    chain_dict = {'Chain ' + letter: {'resids': [], 'angle' : []} for letter in chain_labels}

    print("HERE")
    print(num_prots, chain_length)

    chain_dict['Chain A']['resids'] = [t for t in x if 1 <= t <= chain_length] # 1 to 123

    chain_dict['Chain B']['resids'] = [t for t in x if (chain_length + 1) <= t <= (2 * chain_length)] # 124 to 246

    chain_dict['Chain C']['resids'] = [t for t in x if (1 + (2 * chain_length)) <= t <= (3 * chain_length)] # 247 to 369

    return chain_dict


num_of_proteins = 3
universe = get_universe(sys.argv[1], sys.argv[2])

prot_length = len(universe.select_atoms("protein and name CA"))

chain_dict = read_stride('stride_file.txt', num_of_proteins, int(prot_length / num_of_proteins), False)

print(chain_dict)

# for i in range(3):
#
# sel_resids = ' or resid '.join(map(str, a_ss_resids))
#
# sel = "name CA and (resid " + str(sel_resids) + ")"


#sel = "name CA and resid 124:246"

#print(sel)

# u, ref_princ_axes, princ_axes, dat = run_traj(sel)

princ_axes, data = run_traj(universe, chain_dict)



fig, ax = plt.subplots(len(data), sharex=True, sharey=True)
col = ['#a6cee3','#1f78b4','#b2df8a','#33a02c']

for i, chain in enumerate(data):

    ax[i].scatter(*zip(*data[chain]['angle']), label=str(chain), c=col[i])
    ax[i].legend()
#plt.scatter(*zip(*dat))

# Add titles, legend, and save
#ax.set(title=" ", xlabel="Time (ns)", ylabel="Tilt Angle (" + r'$^{\circ}$' + ")")

fig.text(0.5, 0.04, 'Time (ns)', ha='center')
fig.text(0.04, 0.5, "Tilt Angle (" + r'$^{\circ}$' + ")", va='center', rotation='vertical')

plt.savefig('protein_tilt_pa_to_z_ALLCA.svg', format='svg')

plt.show()


#u = get_universe(sys.argv[1])#, sys.argv[2])















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
#
#
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
#
# print(sys.argv)
# pdb_name = check_argument(sys.argv)
#
# axes = princ_axes
# print(axes)
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
# # print("coordinates:", axis3)
# # print("eigen value:", eval3)
#
# #     print("\nYou can view principal axes with PyMOL:")
# # print("pymol %s %s" %(pymol_name, pdb_name))
#
#
#
#
#
#
#
#
#
#
#
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