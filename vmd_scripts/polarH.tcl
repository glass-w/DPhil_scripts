#Selects the polar hydrogen atoms consistent with OPLS in gromacs
#Source script in TkConsole in VMD (or in a script if you are making representations using a script)
# source ~/pathtoscript/polarH.tcl
#Use in a representation as e.g. "(noh or polarH) and resid 45 47 40"
#If you want to show polar H on backbone use e.g. "(noh or polarHbb) and resid 45 47 40"

#Maria Musgaard, SBCB, April 2013

set def.polarH {((resname SER and name HG) or (resname THR and name HG1) or (resname CYS and name HG) or (resname TYR and name HH) or (resname TRP and name HE1) or (resname ASN and name HD21 HD22) or (resname GLN and name HE21 HE22) or (resname HIS and name HD1 HE2) or (resname LYS and name HZ1 HZ2 HZ3) or (resname ARG and name HE HH11 HH12 HH21 HH22) or (resname GLU and name HE2) or (resname ASP and name HD2))}
atomselect macro polarH ${def.polarH} 

set def.polarHbb {((resname SER and name HG) or (resname THR and name HG1) or (resname CYS and name HG) or (resname TYR and name HH) or (resname TRP and name HE1) or (resname ASN and name HD21 HD22) or (resname GLN and name HE21 HE22) or (resname HIS and name HD1 HE2) or (resname LYS and name HZ1 HZ2 HZ3) or (resname ARG and name HE HH11 HH12 HH21 HH22) or (name H) or (resname GLU and name HE2) or (resname ASP and name HD2))}
atomselect macro polarHbb ${def.polarHbb}
