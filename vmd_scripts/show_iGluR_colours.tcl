#################################################
#						#
#	Script for colouring iGluR and 		#
# 	surroundings				#
#						#
#	Maria Musgaard				#
#	Feb 2012				#
# 	SBCB, Oxford				#
#						#
#################################################

proc colour_iGluR_TMD { molid resolution } {
	# molid: id for molecule to be "visualised"
	# resolution: resolution for representations
	
	source /biggin/s45/musgaard/SCRIPTS/vmd_colourdefs.tcl
	source /biggin/s45/musgaard/SCRIPTS/iGluR_defs.tcl

	set sys $molid
	set res $resolution
	
	# Protein	
	mol representation "NewCartoon 0.250000 $res 2.9500000"
        mol selection "xraypreM1"
	mol material Opaque
        mol color $yellow
        mol addrep $sys
	mol representation "NewCartoon 0.250000 $res 2.9500000"
        mol selection "xrayM1"
	mol material Opaque
        mol color $orange
        mol addrep $sys
	mol representation "NewCartoon 0.250000 $res 2.9500000"
        mol selection "xrayM2"
	mol material Opaque
        mol color $green
        mol addrep $sys
	mol representation "NewCartoon 0.250000 $res 2.9500000"
        mol selection "xrayM3"
	mol material Opaque
        mol color $mauve
        mol addrep $sys
	mol representation "NewCartoon 0.250000 $res 2.9500000"
        mol selection "TMsimM4"
	mol material Opaque
        mol color $cyan
        mol addrep $sys
	mol representation "NewCartoon 0.250000 $res 2.9500000"
        mol selection "protein and not xraypreM1 and not xrayM1 and not xrayM2 and not xrayM3 and not TMsimM4 and not filter"
	mol material Opaque
        mol color $silver
        mol addrep $sys
	mol representation "NewCartoon 0.250000 $res 2.9500000"
        mol selection "filter"
	mol material Opaque
        mol color $blue
        mol addrep $sys
	
	# Specific residues
	mol representation "Licorice 0.300000 $res $res"
        mol selection "noah and protein and resid 617 and (sidechain or name CA)"
       	mol material EdgyShiny
        mol color Name
        mol addrep $sys

	#mol representation "VDW 1.0 $res"
        #mol selection "name CAL"
        #mol material Glossy
        #mol color $cyan3
        #mol addrep $syst
	


	




}
