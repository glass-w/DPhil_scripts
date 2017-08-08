import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import helanal as hel

mda.start_logging()

u = mda.Universe("confout_pbc_corrected.pdb", "traj_pbc_corrected.xtc")

sel = "resnum 10:39 and name CA"


hel.helanal_trajectory(u, sel)

