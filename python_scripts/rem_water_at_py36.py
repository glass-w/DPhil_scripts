import MDAnalysis as mda
import sys

gro_file = sys.argv[1]

u = mda.Universe(gro_file)

p_ul = []
p_ll = []

sol_to_remove = []

sol_com_list = []

p_particles = u.select_atoms("name P")

p_particles_pos = u.select_atoms("name P").positions

#sol_particles = u.select_atoms("resname SOL")
sol_particles = u.select_atoms("name OW")

#sol_particles_resids_com_dict = dict.fromkeys(u.select_atoms("resname SOL"))
sol_particles_resids_com_dict = dict.fromkeys(u.select_atoms("resname SOL and name OW"))


# Find the average z position of the upper and lower leaflets based on P position of lipid
for i, particle in enumerate(p_particles):

    if p_particles_pos[i][2] >= (0.5 * u.dimensions[2]):

        p_ul.append(p_particles_pos[i][2])

    elif p_particles_pos[i][2] <= (0.5 * u.dimensions[2]):

        p_ll.append(p_particles_pos[i][2])

ul = sum(p_ul) / float(len(p_ul))
ll = sum(p_ll) / float(len(p_ll))



# Get the z position of the centre of mass of each water molecule
i = 0
for water in sol_particles:
    print(i)
    #calcualte the CoM of a water molecule and take the z coordinate
    #water_z = (sol_particles[i:i+3].atoms.center_of_mass())[2]

    water_z = (sol_particles.positions[i])[2]

    #assign the z coordinate to the O and both H's of the water molecule for extraction later
    sol_particles_resids_com_dict[sol_particles[i]] = water_z
    #sol_particles_resids_com_dict[sol_particles[i+1]] = water_z
    #sol_particles_resids_com_dict[sol_particles[i + 2]] = water_z

    #i += 3
    i += 1

    #if i + 3 == len(sol_particles) / 4:
     #   print("Progress: 25 %")

   # elif i + 3 == len(sol_particles) / 2:
    #    print("Progress: 50 %")

  #  elif i + 3 == (3 * len(sol_particles)) / 4:
   #     print("Progress: 75 %")

  #  elif i + 3 == len(sol_particles):
  #      break

   # else:
  #      continue

#print(sol_particles_resids_com_dict)

#extract the resid's of the water molecules who's z coordinate of their CoM is not in the bilayer.
sol_resids_not_in_bilayer = [k for k, v in sol_particles_resids_com_dict.items() if v <= ll or v >= ul]

#select everything but water
new_gro_file = u.select_atoms('protein or resname POPC or name NA or name CL')

#make new selection of protein, lipid, ions and water outside of the bilayer
#for group in sol_resids_not_in_bilayer:
for atom in sol_resids_not_in_bilayer:

    new_gro_file += u.select_atoms("resid " + str(atom.resid))

#write out new selection
new_gro_file.write(str(gro_file.split('.')[0]) + ".gro")
