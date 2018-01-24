
proc colorscale_W2R { {reverse 0} {count 0} } {
  display update off
  set mincolorid [expr [colorinfo num] - 1]
  set maxcolorid [expr [colorinfo max] - 1] 
  set colrange [expr $maxcolorid - $mincolorid]
  set colhalf [expr $colrange / 2]

  for {set i $mincolorid} {$i < $maxcolorid} {incr i} {
    set x [expr ($i - $mincolorid) / double($colrange)]

    # quantize so we only get count unique colors, regardless
    # of the low level colorscale matrix size
    if { $count != 0 } {
      set nx [expr {int($count * $x) / double($count)}] 
      set x $nx
    }

    set r 1.0 ;
    set g [expr 1.0 - $x] ;
    set b [expr 1.0 - $x];

    if { $reverse } {
      color change rgb [expr $mincolorid + ($maxcolorid - $i)] $r $g $b 
    } else { 
      color change rgb $i $r $g $b 
		puts "$count $i $r $g $b"
    }
  }

  display update ui
  display update on
}
