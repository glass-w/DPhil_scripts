# Selects the four voltage sensors in CG rep of Nav1.5 

set def.cgvsd1 {resid 90 to 193 and not name W NA CL}
atomselect macro cgvsd1 ${def.cgvsd1}

set def.cgvsd2 {resid 396 to 510 and not name W NA CL}
atomselect macro cgvsd2 ${def.cgvsd2}

set def.cgvsd3 {resid 650 to 760 and not name W NA CL}
atomselect macro cgvsd3 ${def.cgvsd3}

set def.cgvsd4 {resid 950 to 1075 and not name W NA CL}
atomselect macro cgvsd4 ${def.cgvsd4}

set def.all_cgvsds {((resid 90 to 193) or (resid 396 to 510) or (resid 650 to 760) or (resid 950 to 1075)) and not name W NA CL}
atomselect macro all_cgvsds ${def.all_cgvsds}
