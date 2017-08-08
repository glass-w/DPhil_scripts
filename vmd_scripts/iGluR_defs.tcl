###################################################################
###################################################################
###								###
###	Defines the structural elements for iGluR according	###
### 	to the numbers in the sequence in Nature 462, 2009	###
### 	and in the HELIX and SHEET sections of 3KG2.pdb and 	###
###	aditionally visual analysis of the structure.		###
###								###
###	Maria Musgaard						###
###	SBCB, Oxford						###
###								###
###################################################################
###################################################################

# B (or Be) denotes beta strand elemets and A (or Al) alpha helix elements

# Syntax:
# set iGluR.ATD {(protein and resid >= XX and resid <= YY)}
# atomselect macro ATD ${iGluR.ATD}


##### Chains #####

set iGluR.CHNA {(protein and chain A)}
atomselect macro CHNA ${iGluR.CHNA}

set iGluR.CHNB {(protein and chain B)}
atomselect macro CHNB ${iGluR.CHNB}

set iGluR.CHNC {(protein and chain C)}
atomselect macro CHNC ${iGluR.CHNC}

set iGluR.CHND {(protein and chain D)}
atomselect macro CHND ${iGluR.CHND}


##### Domains #####

set iGluR.xrayATD {(protein and resid >= 10 and resid <= 390)}
atomselect macro xrayATD ${iGluR.xrayATD}
 
set iGluR.xrayLBD {(protein and ((resid >=  391 and resid <= 509) or (resid >= 632 and resid <= 772 )))}
atomselect macro xrayLBD ${iGluR.xrayLBD}

set iGluR.xrayTMD {(protein and ((resid >= 510 and resid <= 631) or (resid >= 773)))}
atomselect macro xrayTMD ${iGluR.xrayTMD}

# Residues in M1, M2, M3 and M4
set iGluR.TMDmem {(protein and ((resid >= 514 and resid <= 544) or (resid >= 568 and resid <= 627) or (resid >= 647 and resid <= 674)))}
atomselect macro TMDmem ${iGluR.TMDmem}

# Residues in M1, M2 and M3
set iGluR.TMDmem_noM4 {(protein and ((resid >= 523 and resid <= 544) or (resid >= 568 and resid <= 585) or (resid >= 597 and resid <= 627)))}
atomselect macro TMDmem_noM4 ${iGluR.TMDmem_noM4}


##### Structural elements #####

set iGluR.xrayBe1 {(protein and resid >= 11 and resid <= 19)}
atomselect macro xrayBe1 ${iGluR.xrayBe1}

set iGluR.xrayAl1 {(protein and resid >= 23 and resid <= 36)}
atomselect macro xrayAl1 ${iGluR.xrayAl1}

set iGluR.xrayBe2 {(protein and resid >= 42 and resid <= 50)}
atomselect macro xrayBe2 ${iGluR.xrayBe2}

set iGluR.xrayAl2 {(protein and resid >= 55 and resid <= 68)}
atomselect macro xrayAl2 ${iGluR.xrayAl2}

set iGluR.xrayBe3 {(protein and resid >= 72 and resid <= 75)}
atomselect macro xrayBe3 ${iGluR.xrayBe3}

set iGluR.xrayAl3 {(protein and resid >= 81 and resid <= 92)}
atomselect macro xrayAl3 ${iGluR.xrayAl3}

set iGluR.xrayBe4 {(protein and resid >= 95 and resid <= 98)}
atomselect macro xrayBe4 ${iGluR.xrayBe4}

set iGluR.xrayBe5 {(protein and resid >= 110 and resid <= 112)}
atomselect macro xrayBe5 ${iGluR.xrayBe5}

set iGluR.xrayAl4 {(protein and resid >= 118 and resid <= 128)}
atomselect macro xrayAl4 ${iGluR.xrayAl4}

set iGluR.xrayBe6 {(protein and resid >= 132 and resid <= 137)}
atomselect macro xrayBe6 ${iGluR.xrayBe6}

set iGluR.xrayAl5 {(protein and resid >= 144 and resid <= 157)}
atomselect macro xrayAl5 ${iGluR.xrayAl5}

set iGluR.xrayBe7 {(protein and resid >= 159 and resid <= 164)}
atomselect macro xrayBe7 ${iGluR.xrayBe7}

set iGluR.xrayAl6 {(protein and resid >= 174 and resid <= 187)}
atomselect macro xrayAl6 ${iGluR.xrayAl6}

set iGluR.xrayBe8 {(protein and resid >= 191 and resid <= 195)}
atomselect macro xrayBe8 ${iGluR.xrayBe8}

set iGluR.xrayAl7 {(protein and resid >= 198 and resid <= 212)}
atomselect macro xrayAl7 ${iGluR.xrayAl7}

set iGluR.xrayBe9 {(protein and resid >= 219 and resid <= 222)}
atomselect macro xrayBe9 ${iGluR.xrayBe9}

set iGluR.xrayEta1 {(protein and resid >= 226 and resid <= 229)}
atomselect macro xrayEta1 ${iGluR.xrayEta1}

set iGluR.xrayEta2 {(protein and resid >= 232 and resid <= 237)}
atomselect macro xrayEta2 ${iGluR.xrayEta2}

set iGluR.xrayBe10 {(protein and resid >= 241 and resid <= 246)}
atomselect macro xrayBe10 ${iGluR.xrayBe10}

set iGluR.xrayAl8 {(protein and resid >= 253 and resid <= 265)}
atomselect macro xrayAl8 ${iGluR.xrayAl8}

set iGluR.xrayAl9 {(protein and resid >= 281 and resid <= 303)}
atomselect macro xrayAl9 ${iGluR.xrayAl9}

set iGluR.xrayAl10 {(protein and resid >= 325 and resid <= 336)}
atomselect macro xrayAl10 ${iGluR.xrayAl10}

set iGluR.xrayBe11 {(protein and resid >= 339 and resid <= 340)}
atomselect macro xrayBe11 ${iGluR.xrayBe11}

set iGluR.xrayBe12 {(protein and resid >= 343 and resid <= 344)}
atomselect macro xrayBe12 ${iGluR.xrayBe12}

set iGluR.xrayBe13 {(protein and resid >= 357 and resid <= 363)}
atomselect macro xrayBe13 ${iGluR.xrayBe13}

set iGluR.xrayBe14 {(protein and resid >= 368 and resid <= 371)}
atomselect macro xrayBe14 ${iGluR.xrayBe14}

set iGluR.xrayBe15 {(protein and resid >= 373 and resid <= 375)}
atomselect macro xrayBe15 ${iGluR.xrayBe15}

set iGluR.xrayB1 {(protein and resid >= 395 and resid <= 399)}
atomselect macro xrayB1 ${iGluR.xrayB1}

set iGluR.xrayaA1 {(protein and resid >= 411 and resid <= 415)}
atomselect macro xrayaA1 ${iGluR.xrayaA1}

set iGluR.xrayaA2 {(protein and resid >= 417 and resid <= 419)}
atomselect macro xrayaA2 ${iGluR.xrayaA2}

set iGluR.xrayaB {(protein and resid >= 424 and resid <= 437)}
atomselect macro xrayaB ${iGluR.xrayaB}

set iGluR.xrayB2 {(protein and resid >= 440 and resid <= 444)}
atomselect macro xrayB2 ${iGluR.xrayB2}

set iGluR.xrayB3 {(protein and resid >= 451 and resid <= 454)}
atomselect macro xrayB3 ${iGluR.xrayB3}

set iGluR.xrayB4 {(protein and resid >= 458 and resid <= 460)}
atomselect macro xrayB4 ${iGluR.xrayB4}

set iGluR.xrayaC {(protein and resid >= 462 and resid <= 469)}
atomselect macro xrayaC ${iGluR.xrayaC}

set iGluR.xrayB5 {(protein and resid >= 474 and resid <= 475)}
atomselect macro xrayB5 ${iGluR.xrayB5}

set iGluR.xrayaD {(protein and resid >= 483 and resid <= 488)}
atomselect macro xrayaD ${iGluR.xrayaD}

set iGluR.xrayB6 {(protein and resid >= 490 and resid <= 491)}
atomselect macro xrayB6 ${iGluR.xrayB6}

set iGluR.xrayB7 {(protein and resid >= 496 and resid <= 498)}
atomselect macro xrayB7 ${iGluR.xrayB7}

set iGluR.xrayB8 {(protein and resid >= 500 and resid <= 505)}
atomselect macro xrayB8 ${iGluR.xrayB8}

set iGluR.xraypreM1 {(protein and resid >= 514 and resid <= 521)}
atomselect macro xraypreM1 ${iGluR.xraypreM1}

set iGluR.xrayM1 {(protein and resid >= 523 and resid <= 544)}
atomselect macro xrayM1 ${iGluR.xrayM1}

set iGluR.xrayM2 {(protein and resid >= 568 and resid <= 585)}
atomselect macro xrayM2 ${iGluR.xrayM2}

set iGluR.xrayM3 {(protein and resid >= 597 and resid <= 627)}
atomselect macro xrayM3 ${iGluR.xrayM3}

set iGluR.xrayaE {(protein and resid >= 636 and resid <= 642)}
atomselect macro xrayaE ${iGluR.xrayaE}

set iGluR.xrayB9 {(protein and resid >= 646 and resid <= 648)}
atomselect macro xrayB9 ${iGluR.xrayB9}

set iGluR.xrayaF {(protein and resid >= 654 and resid <= 662)}
atomselect macro xrayaF ${iGluR.xrayaF}

set iGluR.xrayaG {(protein and resid >= 665 and resid <= 677)}
atomselect macro xrayaG ${iGluR.xrayaG}

set iGluR.xrayaH {(protein and resid >= 686 and resid <= 696)}
atomselect macro xrayaH ${iGluR.xrayaH}

set iGluR.xrayB10 {(protein and resid >= 700 and resid <= 705)}
atomselect macro xrayB10 ${iGluR.xrayB10}

set iGluR.xrayaI {(protein and resid >= 706 and resid <= 715)}
atomselect macro xrayaI ${iGluR.xrayaI}

set iGluR.xrayB11 {(protein and resid >= 720 and resid <= 724)}
atomselect macro xrayB11 ${iGluR.xrayB11}

set iGluR.xrayB12 {(protein and resid >= 727 and resid <= 731)}
atomselect macro xrayB12 ${iGluR.xrayB12}

set iGluR.xrayB13 {(protein and resid >= 735 and resid <= 737)}
atomselect macro xrayB13 ${iGluR.xrayB13}

set iGluR.xrayaJ {(protein and resid >= 742 and resid <= 756)}
atomselect macro xrayaJ ${iGluR.xrayaJ}

set iGluR.xrayaK {(protein and resid >= 758 and resid <= 769)}
atomselect macro xrayaK ${iGluR.xrayaK}

set iGluR.xraypreM4 {(protein and resid >= 774 and resid <= 781)}
atomselect macro xraypreM4 ${iGluR.xraypreM4}

set iGluR.xrayM4 {(protein and resid >= 790 and resid <= 817)}
atomselect macro xrayM4 ${iGluR.xrayM4}

set iGluR.TMsimM4 {(protein and resid >= 647 and resid <= 674)}
atomselect macro TMsimM4 ${iGluR.TMsimM4}

set iGluR.xraySS {(protein and resid 63 315 718 773)}
atomselect macro xraySS ${iGluR.xraySS}

set iGluR.xraydomS1 {(protein and resid >= 391 and resid <=509)}
atomselect macro xraydomS1 ${iGluR.xraydomS1}

set iGluR.xraydomS2 {(protein and resid >= 632 and resid <=772)}
atomselect macro xraydomS2 ${iGluR.xraydomS2}

set iGluR.filter {(protein and resid >= 586 and resid <=592)}
atomselect macro filter ${iGluR.filter}

