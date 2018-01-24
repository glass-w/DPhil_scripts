# selects the TM region of beta3 monomer in CG rep

set def.cg_beta3_tm {resid 136 to 156 and not name W NA CL}
atomselect macro cg_beta3_tm ${def.cg_beta3_tm}
