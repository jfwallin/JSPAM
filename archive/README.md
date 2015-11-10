Archive
=====
This directory allows you to compile the simulation code in a form that is almost identical to what was run in the applet used in Galaxy Zoo: Mergers.  The primary changes were to remove the dependence on the web server and to allow for simpler execution.

Compiling
---------
To compile, simply use make or ant.

Execution
---------
After the code has been compiled, run the runall.sh script.  It will create an output directory which will then have a single directory for each of the 62 targets.  The simulation can be plotted in gnuplot by executing gnuplot < plot from each directory.
