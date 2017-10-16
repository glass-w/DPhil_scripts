#***************************
#show version and usage info
#***************************
proc vmd_parser_info {} {
	set version "0.1.0"
	puts "\[ ABOUT \]"	
	puts "version: $version"
	puts "author: Jean Helie (jean.helie@bioch.ox.ac.uk)"
	puts "github: https://github.com/jhelie/vmd_parser"
	puts ""
	puts "\[ USAGE \]"
	puts "Type 'function user_field file.txt', where 'function' can be: "
	puts " -> set_cluster_prot"
	puts " -> set_cluster_lip"
	puts " -> set_order_param"
	puts " -> set_thickness"
	puts " -> set_leaflet"
	puts "and where 'user_field' can be user, user2, user3 or user4 (in the case of set_leaflet no need to specify a user field, the info is stored in the beta field)."
	puts ""
}

#***************************************************************************
#load protein clustering info outputed by the MDAnalysis script cluster_prot
#***************************************************************************
proc set_cluster_prot {txtfilename userfield} {
	
	#check userfield option
	#----------------------
	set userfield_values {user user2 user3 user4}
	if {$userfield in $userfield_values} {
		puts "data stored in : $userfield"
	} else {
		puts "Error incorrect option for the user field, see vmd_parser_info"
		return
	}
	
	#read file contents
	#------------------
	set file_handle [open $txtfilename]
	set file_content [read $file_handle]
	set lines [split $file_content "\n"]
	set lines [lreplace $lines [expr [llength $lines] - 1] [expr [llength $lines] - 1]]
		
	#store header line and remove it from list
	#-----------------------------------------
	set prot_sele_string [split [lindex $lines 0] .]
	set lines [lreplace $lines 0 0]
	set prot_nb [llength $prot_sele_string]
	set prot_indexes {0}
	for {set i 1} {$i<$prot_nb} {incr i} {
		set prot_indexes [lappend prot_indexes $i]
	}

	#store min and max sizes line and remove it from list
	#----------------------------------------------------
	set sizes [split [lindex $lines 0] .]
	set lines [lreplace $lines 0 0]
	set g_min [lindex $sizes 0]
	set g_max [lindex $sizes 1]

	#check coherence between length of xtc and of data file
	#------------------------------------------------------
	set nb_lines [llength $lines]
	set nb_frames [molinfo top get numframes]
	if {$nb_frames!=1 && $nb_lines != [expr $nb_frames -1]} {			#the -1 is because the gro/pdb file is counted as a frame
		puts "Warning: different nb of frames ([expr $nb_frames -1]) and nb of lines ($nb_lines)."
		puts "Only the frames referenced in the annotation file will be annotated."
	} 

	#process each frame
	#------------------
	if {$userfield == "user"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line .]
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant cluster info for each protein
			foreach p_index $prot_indexes { 
				set tmp_sele [atomselect top [lindex $prot_sele_string $p_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user [lindex $tmp_line [expr $p_index +1] ]
				$tmp_sele delete
			}
		}
	} elseif {$userfield == "user2"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line .]
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant cluster info for each protein
			foreach p_index $prot_indexes { 
				set tmp_sele [atomselect top [lindex $prot_sele_string $p_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user2 [lindex $tmp_line [expr $p_index +1] ]
				$tmp_sele delete
			}
		}
	} elseif {$userfield == "user3"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line .]
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant cluster info for each protein
			foreach p_index $prot_indexes { 
				set tmp_sele [atomselect top [lindex $prot_sele_string $p_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user3 [lindex $tmp_line [expr $p_index +1] ]
				$tmp_sele delete
			}
		}
	} elseif {$userfield == "user4"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line .]
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant cluster info for each protein
			foreach p_index $prot_indexes { 
				set tmp_sele [atomselect top [lindex $prot_sele_string $p_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user4 [lindex $tmp_line [expr $p_index +1] ]
				$tmp_sele delete
			}
		}	
	}
	
	
	#remind user of coloring options
	#-------------------------------
	puts "Done. To color proteins by cluster size status:"
	puts " 1.select Draw Style > Coloring Method > Trajectory > User > $userfield"
	puts " 2.set Trajectory > Color Scale Data Range: min=$g_min, max=$g_max and click Set."
	puts " 3.set Trajectory > Update Color Every Frame"
	puts ""
	puts "NB: if you use custom colorscale from matlab_colorscale.tcl you need want to change \"'\$i < \$maxcolorid'\" to \"'\$i <= \$maxcolorid'\" in the file in order to reproduce the correct maximum colour."
	puts ""
}

proc wsplit {str sep} {
  split [string map [list $sep \0] $str] \0
  }
  
#***************************************************************************
#load protein clustering info outputed by the meso model
#***************************************************************************
proc set_cluster_prot_meso {txtfilename userfield nprot} {
	
	#check userfield option
	#----------------------
	set userfield_values {user user2 user3 user4}
	if {$userfield in $userfield_values} {
		puts "data stored in : $userfield"
	} else {
		puts "Error incorrect option for the user field, see vmd_parser_info"
		return
	}
	
	#read file contents
	#------------------
	set file_handle [open $txtfilename]
	set file_content [read $file_handle]
	set lines [split $file_content "\n"]
	set lines [lreplace $lines [expr [llength $lines] - 1] [expr [llength $lines] - 1]]
		
	#store header line and remove it from list
	#-----------------------------------------
	#set prot_sele_string [split [lindex $lines 0] .]
	#set lines [lreplace $lines 0 0]
	#set prot_nb [llength $prot_sele_string]
	#set prot_indexes {0}
	#for {set i 1} {$i<$prot_nb} {incr i} {
	#	set prot_indexes [lappend prot_indexes $i]
	#}

	#store min and max sizes line and remove it from list
	#----------------------------------------------------
	#set sizes [split [lindex $lines 0] .]
	#set lines [lreplace $lines 0 0]
	#set g_min [lindex $sizes 0]
	#set g_max [lindex $sizes 1]

	#check coherence between length of xtc and of data file
	#------------------------------------------------------
	set nb_lines [llength $lines]
	set nb_frames [molinfo top get numframes]
	if {$nb_frames!=1 && $nb_lines != [expr $nb_frames -1]} {			#the -1 is because the gro/pdb file is counted as a frame
		puts "Warning: different nb of frames ([expr $nb_frames -1]) and nb of lines ($nb_lines)."
		puts "Only the frames referenced in the annotation file will be annotated."
	} 

	#process each frame
	#------------------
	if {$userfield == "user"} {
	    set tmp_frame 0
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line0 [string trim $line "\[\]"]
			set tmp_line [wsplit $tmp_line0 ", "]
			set tmp_frame [ expr $tmp_frame + 1 ]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant cluster info for each protein
			for {set p_index 1} {$p_index<$nprot} {incr p_index} {
				set tmp_sele [atomselect top "serial $p_index"]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user [lindex $tmp_line $p_index ]
				$tmp_sele delete
			}
		}
	} elseif {$userfield == "user2"} {
	    set tmp_frame 0
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line0 [string trim $line "\[\]"]
			set tmp_line [wsplit $tmp_line0 ", "]
			set tmp_frame [ expr $tmp_frame + 1 ]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant cluster info for each protein
			foreach p_index $prot_indexes { 
				set tmp_sele [atomselect top "serial $i"]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user2 [lindex $tmp_line $p_index ]
				$tmp_sele delete
			}
		}
	} elseif {$userfield == "user3"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line .]
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant cluster info for each protein
			foreach p_index $prot_indexes { 
				set tmp_sele [atomselect top [lindex $prot_sele_string $p_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user3 [lindex $tmp_line [expr $p_index +1] ]
				$tmp_sele delete
			}
		}
	} elseif {$userfield == "user4"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line .]
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant cluster info for each protein
			foreach p_index $prot_indexes { 
				set tmp_sele [atomselect top [lindex $prot_sele_string $p_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user4 [lindex $tmp_line [expr $p_index +1] ]
				$tmp_sele delete
			}
		}	
	}
	
	
	#remind user of coloring options
	#-------------------------------
	puts "Done. To color proteins by cluster size status:"
	puts " 1.select Draw Style > Coloring Method > Trajectory > User > $userfield"
	puts " 2.set Trajectory > Color Scale Data Range: min=$g_min, max=$g_max and click Set."
	puts " 3.set Trajectory > Update Color Every Frame"
	puts ""
	puts "NB: if you use custom colorscale from matlab_colorscale.tcl you need want to change \"'\$i < \$maxcolorid'\" to \"'\$i <= \$maxcolorid'\" in the file in order to reproduce the correct maximum colour."
	puts ""
}
#************************************************************************
#load lipid clustering info outputed by the MDAnalysis script cluster_lip
#************************************************************************
proc set_cluster_lip {txtfilename userfield} {
	
	#check userfield option
	#----------------------
	set userfield_values {user user2 user3 user4}
	if {$userfield in $userfield_values} {
		puts "data stored in : $userfield"
	} else {
		puts "Error incorrect option for the user field, see vmd_parser_info"
		return
	}

	#read file contents
	#------------------
	set file_handle [open $txtfilename]
	set file_content [read $file_handle]
	set lines [split $file_content "\n"]
	set lines [lreplace $lines [expr [llength $lines] - 1] [expr [llength $lines] - 1]]
		
	#store header lines and remove them from list
	#--------------------------------------------
	set lip_sele_string [split [lindex $lines 0] .]
	set lines [lreplace $lines 0 0]
	set lip_nb [llength $lip_sele_string]
	set lip_indexes {0}
	for {set i 1} {$i<$lip_nb} {incr i} {
		set lip_indexes [lappend lip_indexes $i]
	}

	#store min and max thickness line and remove it from list
	#--------------------------------------------------------
	set range [split [lindex $lines 0] ";"]
	set lines [lreplace $lines 0 0]
	set c_min [lindex $range 0]
	set c_max [lindex $range 1]
	
	#check coherence between length of xtc and of data file
	#------------------------------------------------------
	set nb_lines [llength $lines]
	set nb_frames [molinfo top get numframes]
	if {$nb_frames!=1 && $nb_lines != [expr $nb_frames -1]} { 			#the -1 is because the gro/pdb file is counted as a frame
		puts "Warning: different nb of frames ([expr $nb_frames -1]) and nb of lines ($nb_lines)."
		puts "Only the frames referenced in the annotation file will be annotated."
	} 

	#process each frame
	#------------------
	if {$userfield == "user"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line ";"]			
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant order param info for each lipid
			foreach r_index $lip_indexes { 
				set tmp_sele [atomselect top [lindex $lip_sele_string $r_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user [lindex $tmp_line [expr $r_index +1] ]
				$tmp_sele delete
			}
		}
	} elseif {$userfield == "user2"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line ";"]			
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant order param info for each lipid
			foreach r_index $lip_indexes { 
				set tmp_sele [atomselect top [lindex $lip_sele_string $r_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user2 [lindex $tmp_line [expr $r_index +1] ]
				$tmp_sele delete
			}
		}
	} elseif {$userfield == "user3"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line ";"]			
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant order param info for each lipid
			foreach r_index $lip_indexes { 
				set tmp_sele [atomselect top [lindex $lip_sele_string $r_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user3 [lindex $tmp_line [expr $r_index +1] ]
				$tmp_sele delete
			}
		}
	} elseif {$userfield == "user4"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line ";"]			
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant order param info for each lipid
			foreach r_index $lip_indexes { 
				set tmp_sele [atomselect top [lindex $lip_sele_string $r_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user4 [lindex $tmp_line [expr $r_index +1] ]
				$tmp_sele delete
			}
		}
	}

	#remind user of coloring options
	#-------------------------------
	puts "Done. To color lipids by cluster size status:"
	puts " 1.select Draw Style > Coloring Method > Trajectory > User > $userfield"
	puts " 2.set Trajectory > Color Scale Data Range: min=$c_min, max=$c_max and click Set."
	puts " 3.set Trajectory > Update Color Every Frame"
	puts ""
	puts "NB: if you use custom colorscale from matlab_colorscale.tcl you need want to change \"'\$i < \$maxcolorid'\" to \"'\$i <= \$maxcolorid'\" in the file in order to reproduce the correct maximum colour."
	puts ""
}

#*******************************************************************************
#load bilayer order parameter info outputed by the MDAnalysis script order_param
#*******************************************************************************
proc set_order_param {txtfilename userfield} {

	#check userfield option
	#----------------------
	set userfield_values {user user2 user3 user4}
	if {$userfield in $userfield_values} {
		puts "data stored in : $userfield"
	} else {
		puts "Error incorrect option for the user field, see vmd_parser_info"
		return
	}
	
	#read file contents
	#------------------
	set file_handle [open $txtfilename]
	set file_content [read $file_handle]
	set lines [split $file_content "\n"]
	set lines [lreplace $lines [expr [llength $lines] - 1] [expr [llength $lines] - 1]]
		
	#store header lines and remove them from list
	#--------------------------------------------
	set lip_sele_string [split [lindex $lines 0] .]
	set lines [lreplace $lines 0 0]
	set lip_nb [llength $lip_sele_string]
	set lip_indexes {0}
	for {set i 1} {$i<$lip_nb} {incr i} {
		set lip_indexes [lappend lip_indexes $i]
	}

	#check coherence between length of xtc and of data file
	#------------------------------------------------------
	set nb_lines [llength $lines]
	set nb_frames [molinfo top get numframes]
	if {$nb_frames!=1 && $nb_lines != [expr $nb_frames -1]} { 			#the -1 is because the gro/pdb file is counted as a frame
		puts "Warning: different nb of frames ([expr $nb_frames -1]) and nb of lines ($nb_lines)."
		puts "Only the frames referenced in the annotation file will be annotated."
	} 

	#process each frame
	#------------------
	if {$userfield == "user"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line ";"]			
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant order param info for each lipid
			foreach r_index $lip_indexes { 
				set tmp_sele [atomselect top [lindex $lip_sele_string $r_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user [lindex $tmp_line [expr $r_index +1] ]
				$tmp_sele delete
			}
		}
	} elseif {$userfield == "user2"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line ";"]			
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant order param info for each lipid
			foreach r_index $lip_indexes { 
				set tmp_sele [atomselect top [lindex $lip_sele_string $r_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user2 [lindex $tmp_line [expr $r_index +1] ]
				$tmp_sele delete
			}
		}
	} elseif {$userfield == "user3"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line ";"]			
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant order param info for each lipid
			foreach r_index $lip_indexes { 
				set tmp_sele [atomselect top [lindex $lip_sele_string $r_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user3 [lindex $tmp_line [expr $r_index +1] ]
				$tmp_sele delete
			}
		}
	} elseif {$userfield == "user4"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line ";"]			
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant order param info for each lipid
			foreach r_index $lip_indexes { 
				set tmp_sele [atomselect top [lindex $lip_sele_string $r_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user4 [lindex $tmp_line [expr $r_index +1] ]
				$tmp_sele delete
			}
		}
	}
		
	#remind user of coloring options
	#-------------------------------
	puts "Done. To color lipids based on their tail order parameter:"
	puts " 1.select Draw Style > Coloring Method > Trajectory > User > $userfield"
	puts " 2.set Trajectory > Color Scale Data Range: min=-0.5, max=1 and click Set."
	puts " 3.set Trajectory > Update Color Every Frame"
	puts ""
	puts "NB: if you use custom colorscale from matlab_colorscale.tcl you need want to change \"'\$i < \$maxcolorid'\" to \"'\$i <= \$maxcolorid'\" in the file in order to reproduce the correct maximum colour."
	puts ""
}

#***********************************************************************
#load bilayer thickness info outputed by the MDAnalysis script thickness
#***********************************************************************
proc set_thickness {txtfilename userfield} {
	#check userfield option
	#----------------------
	set userfield_values {user user2 user3 user4}
	if {$userfield in $userfield_values} {
		puts "data stored in : $userfield"
	} else {
		puts "Error incorrect option for the user field, see vmd_parser_info"
		return
	}
	
	#read file contents
	#------------------
	set file_handle [open $txtfilename]
	set file_content [read $file_handle]
	set lines [split $file_content "\n"]
	set lines [lreplace $lines [expr [llength $lines] - 1] [expr [llength $lines] - 1]]
		
	#store header lines and remove them from list
	#--------------------------------------------
	set lip_sele_string [split [lindex $lines 0] .]
	set lines [lreplace $lines 0 0]
	set lip_nb [llength $lip_sele_string]
	set lip_indexes {0}
	for {set i 1} {$i<$lip_nb} {incr i} {
		set lip_indexes [lappend lip_indexes $i]
	}

	#store min and max thickness line and remove it from list
	#--------------------------------------------------------
	set range [split [lindex $lines 0] ";"]
	set lines [lreplace $lines 0 0]
	set h_min [lindex $range 0]
	set h_max [lindex $range 1]

	#check coherence between length of xtc and of data file
	#------------------------------------------------------
	set nb_lines [llength $lines]
	set nb_frames [molinfo top get numframes]
	if {$nb_frames!=1 && $nb_lines != [expr $nb_frames -1]} { 			#the -1 is because the gro/pdb file is counted as a frame
		puts "Warning: different nb of frames ([expr $nb_frames -1]) and nb of lines ($nb_lines)."
		puts "Only the frames referenced in the annotation file will be annotated."
	} 

	#process each frame
	#------------------
	if {$userfield == "user"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line ";"]			
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant order param info for each lipid
			foreach r_index $lip_indexes { 
				set tmp_sele [atomselect top [lindex $lip_sele_string $r_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user [lindex $tmp_line [expr $r_index +1] ]
				$tmp_sele delete
			}
		}
	} elseif {$userfield == "user2"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line ";"]			
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant order param info for each lipid
			foreach r_index $lip_indexes { 
				set tmp_sele [atomselect top [lindex $lip_sele_string $r_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user2 [lindex $tmp_line [expr $r_index +1] ]
				$tmp_sele delete
			}
		}
	} elseif {$userfield == "user3"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line ";"]			
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant order param info for each lipid
			foreach r_index $lip_indexes { 
				set tmp_sele [atomselect top [lindex $lip_sele_string $r_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user3 [lindex $tmp_line [expr $r_index +1] ]
				$tmp_sele delete
			}
		}
	} elseif {$userfield == "user4"} {
		foreach line $lines { 
			#retrieve current line and the frame it references
			set tmp_line [split $line ";"]			
			set tmp_frame [lindex $tmp_line 0]
			
			#display
			puts "processing frame $tmp_frame..."
			
			#set relevant order param info for each lipid
			foreach r_index $lip_indexes { 
				set tmp_sele [atomselect top [lindex $lip_sele_string $r_index]]
				$tmp_sele frame $tmp_frame
				$tmp_sele set user4 [lindex $tmp_line [expr $r_index +1] ]
				$tmp_sele delete
			}
		}
	}
		
	#remind user of coloring options
	#-------------------------------
	puts "Done. To color lipids based on the bilayer thickness:"
	puts " 1.select Draw Style > Coloring Method > Trajectory > User > $userfield"
	puts " 2.set Trajectory > Color Scale Data Range: min=$h_min, max=$h_max and click Set."
	puts " 3.set Trajectory > Update Color Every Frame"
	puts ""
	puts "NB: if you use custom colorscale from matlab_colorscale.tcl you need want to change \"'\$i < \$maxcolorid'\" to \"'\$i <= \$maxcolorid'\" in the file in order to reproduce the correct maximum colour."
	puts ""
}

#***********************************************************************
#load leaflet info outputed by the MDAnalysis script leaflet_annotator
#***********************************************************************
proc set_leaflet {txtfilename} {
	
	#read file contents
	#------------------
	set file_handle [open $txtfilename]
	set file_content [read $file_handle]
	set lines [split $file_content "\n"]
	set lines [lreplace $lines [expr [llength $lines] - 1] [expr [llength $lines] - 1]]
		
	#store header lines and remove them from list
	#--------------------------------------------
	set lip_sele_string [split [lindex $lines 0] .]
	set lines [lreplace $lines 0 0]
	set lip_nb [llength $lip_sele_string]
	set lip_indexes {0}
	for {set i 1} {$i<$lip_nb} {incr i} {
		set lip_indexes [lappend lip_indexes $i]
	}

	#set beta factor (leaflet aren't changing through xtc, no need to use the user fields..."
	#---------------
	foreach line $lines { 
		#retrieve current line and the frame it references
		set tmp_line [split $line .]			
		set tmp_frame [lindex $tmp_line 0]
					
		#set relevant order param info for each lipid
		foreach r_index $lip_indexes { 
			set tmp_sele [atomselect top [lindex $lip_sele_string $r_index]]
			$tmp_sele set beta [lindex $tmp_line [expr $r_index +1] ]
			$tmp_sele delete
		}
	}
		
	#remind user of coloring options
	#-------------------------------
	puts "Done. To color lipids by leaflet status:"
	puts " select Draw Style > Coloring Method > Beta"
	puts ""
}

