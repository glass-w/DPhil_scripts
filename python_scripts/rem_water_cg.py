import MDAnalysis as mda

u = mda.Universe("b3_6x6_system_red_z_solvated.gro")

po4_ul = []
po4_ll = []

w_to_remove = []

po4_particles = u.select_atoms("name PO4")

po4_particles_pos = u.select_atoms("name PO4").positions

w_particles = u.select_atoms("name W")

w_particles_pos = u.select_atoms("name W").positions

for particle in range(len(po4_particles)):

    if po4_particles_pos[particle][2] >= (0.5 * u.dimensions[2]):

        po4_ul.append(po4_particles_pos[particle][2])

    elif po4_particles_pos[particle][2] < (0.5 * u.dimensions[2]):

        po4_ll.append(po4_particles_pos[particle][2])

ul = sum(po4_ul) / float(len(po4_ul))
ll = sum(po4_ll) / float(len(po4_ll))

w_above_ul = u.select_atoms("name W and prop z > " + str(ul))

w_below_ul = u.select_atoms("name W and prop z < " + str(ll))

test = w_above_ul + w_below_ul


new_gro_file = u.select_atoms("all and not (name W and (prop z < " + str(ul) + " and prop z > " + str(ll) + "))")

#below = u.select_atoms("all and not (name W and prop z > " + str(ll) + ")")

#w_removed = above + below

#print (int(len(u.select_atoms("resname POPC"))) / 12)

#w_removed.write("b3_6x6_system_red_z_solvated.gro")

new_gro_file.write("b3_6x6_system_red_z_solvated.gro")
