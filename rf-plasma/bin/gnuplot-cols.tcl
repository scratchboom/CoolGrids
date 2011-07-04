#!/usr/bin/tclsh
# Auxilary procedures
proc vval {v} {return $v}

set slash "/"
proc base_name {path} {
	global slash
	set f [string last $slash $path]
	if {$f == -1} {
		set f 0
	} else {
		set f [expr $f + 1]
	}
	set l [string last "." $path]
	if {$l == -1} {
		set l [expr [string length $path] - 1]
	} else {
		set l [expr $l - 1]
	}
	return [string range $path $f $l]
}

proc full_name {path} {
	global slash
	set f [string last $slash $path]
	if {$f == -1} {
		set f 0
	} else {
		set f [expr $f + 1]
	}
	set l [expr [string length $path] - 1]
	return [string range $path $f $l]
}

proc parent_dir {path} {
	global slash
	set l [string last $slash $path]
	if {$l == -1} {
		return ""
	} else {
		set l [expr $l - 1]
	}
	return [string range $path 0 $l]
}

proc dir_attach {path child} {
	global slash
	if {[string length $path] > 0} {
		return "[vval $path]$slash[vval $child]"
	} else {
		return $child
	}
}

proc generate_plot_arguments {infile_l} {
        global MODE
	set result ""
	set i 0
	set filename {}
	set cols {}
	set title {}
	foreach elem $infile_l {
		if {[expr $i % 3] == 0} {
			set filename $elem
		}
		if {[expr $i % 3] == 1} {
			set cols $elem
		}
		if {[expr $i % 3] == 2} {
			if {$elem == "-"} {
				set title [base_name $filename]-$cols
			} else {
				set title $elem
			}
			set result "$result, '$filename' using $cols with $MODE title \"$title\""
		}
		incr i
	}
	return [string range $result 2 [expr [string length $result] - 1]]
}

proc generate_plot_title {title} {
	if {[string length $title] != 0} {
		return "set title \"$title\""
	} else {
		return ""
	}
}

# Global variables
set help_text "\
Usage: tclsh [full_name $argv0] INPUT_FILES \[OUTPUT_FILE\] \[EXTRA_COMMANDS\] \[TITLE\] \[MODE\]
      \[XLABEL\] \[YLABEL\] \[XSIZE,YSIZE\] \[GNUPLOT_OUTPUT_FILE\]

 ,where INPUT_FILES = \"FILENAME1 XCOL1:YCOL1 TITLE1
                    \[FILENAME2 XCOL2:YCOL2 TITLE2\] ... \"
Generates gnuplot script and executes it
"

# Command-line arguments parsing
if {$argc < 1} {
	puts $help_text
	exit
}

# Requied parameters
set INPUT_FILES [lindex $argv 0]

# Defaults for optional parameters
if {$argc == 1} {
	if {[llength $INPUT_FILES] != 1} {
		puts "ERROR: The second argument should be provided."
		exit
	}
	set OUTPUT_FILE [dir_attach [parent_dir $INPUT_FILES] [base_name $INPUT_FILES].png]
}

set TITLE ""
set MODE "lines"
set XLABEL ""
set YLABEL ""
set SIZE "1400,1050"
set EXTRA_COMMANDS ""
set GNUPLOT_OUTPUT_FILE "temp.[exec perl -e {print $$}].gp"

set delete_gp true

# Optional parameters
if {$argc > 1} {set OUTPUT_FILE [lindex $argv 1]}
if {$argc > 2} {if {[lindex $argv 2] != "-"} {set EXTRA_COMMANDS [lindex $argv 2]} }
if {$argc > 3} {set TITLE [lindex $argv 3]}
if {$argc > 4} {set MODE [lindex $argv 4]}
if {$argc > 5} {set XLABEL [lindex $argv 5]}
if {$argc > 6} {set YLABEL [lindex $argv 6]}
if {$argc > 7} {set SIZE [lindex $argv 7]}
if {$argc > 8} {set GNUPLOT_OUTPUT_FILE [lindex $argv 8];set delete_gp false}

# Generating gnuplot script (was developed under gnuplot version 4.2.0)
set gpf [open $GNUPLOT_OUTPUT_FILE w]
puts $gpf "

# Parametes
[generate_plot_title $TITLE]
#set autotitle
set terminal png size $SIZE
set output \"$OUTPUT_FILE\"
set xlabel \"$XLABEL\"
set ylabel \"$YLABEL\"
#set datafile commentschars \"\\\"\"

# Extra commands
$EXTRA_COMMANDS

# Plotting
plot [generate_plot_arguments $INPUT_FILES];
"
close $gpf

# Running script
exec gnuplot $GNUPLOT_OUTPUT_FILE

# Remove script if needed
if {$delete_gp} {
	exec rm -f $GNUPLOT_OUTPUT_FILE
}

