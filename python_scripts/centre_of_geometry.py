import MDAnalysis as mda
import numpy as np
from matplotlib import pyplot as plt
import argparse

def cog(gro_file, traj_file):

     u = mda.Universe(gro_file, traj_file)

     # define regions of which you want to calculate the centre of geometry
     # ig_a = u.select_atoms("resid 1:123")
     # ig_b = u.select_atoms("resid 167:290")
     # ig_c = u.select_atoms("resid 333:457")

     ig_a = u.select_atoms("resid 1:123")
     ig_b = u.select_atoms("resid 124:246")
     ig_c = u.select_atoms("resid 247:369")

     r1_store = []
     r2_store = []
     r3_store = []
     time = []

     for frame in range(0, len(u.trajectory), 10):
          print frame
          u.trajectory[frame]
          time.append((u.trajectory.time)/1000)
          print ((u.trajectory.time)/1000)

          ig_a_com = ig_a.select_atoms("name CA").center_of_mass()
          ig_b_com = ig_b.select_atoms("name CA").center_of_mass()
          ig_c_com = ig_c.select_atoms("name CA").center_of_mass()

          r1 = np.linalg.norm(ig_a_com - ig_c_com)
          r1_store.append(r1)

          r2 = np.linalg.norm(ig_a_com - ig_b_com)
          r2_store.append(r2)

          r3 = np.linalg.norm(ig_b_com - ig_c_com)
          r3_store.append(r3)

     ax = plt.subplot(111)
     ax.set_xlabel("Time (ns)")
     ax.set_ylabel("Distance ($\AA$)")

     ax.plot(time, r1_store, 'b', label='ig_A - ig_C', alpha=0.4)
     ax.plot(time, r2_store, 'r', label='ig_A - ig_B', alpha=0.4)
     ax.plot(time, r3_store, 'g', label='ig_B - ig_C', alpha=0.4)
     ax.legend()
     plt.title("Distances Between Subunits EC Domains (water box)")
     plt.savefig('CoG_plot_b3TrimerWB_r2.svg', dpi=300)
     plt.show()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculates the centre of geometry between three defined protein regions")
    parser.add_argument("-c", dest="gro_file", type=str, help='the coordinate file [.gro]')
    parser.add_argument("-f", dest="xtc_file", type=str,
                        help='a corrected trajectory file, pbc artifacts removed and the protein centered')
    options = parser.parse_args()

    cog(options.gro_file, options.xtc_file)

