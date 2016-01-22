########################################################################
#                gnuplot script exactMovie.gnu
#
# A script for creating an animation of the data files created by 
# Assign1.cpp. This script must be executed in the directory containiing
# the output data files uExactXXX.dat where XXX is the output number. 
#
# If the maximal output number is N, then to invoke the script use the 
# sequence of commands
#
# > output_count = 0 
# > output_max   = N
# > load "exactMovie.gnu"
#
# where N is replaced with the maximal output number e.g. 50.
#
# The variable output_count has to be reset in order to re-run
# the script. 
#
# To view the animated gif in uExact.gif, open the file 
# with your web browser (Firefox works well) or use QuickTime. 
# With the browser, hitting the ESC key will stop the animation. 
#
# Created for Math 269B 
# Creation Data Jan. 2, 2014
# Version: Jan. 6, 2016 
# Author Chris Anderson

########################################################################

#
# For the first iteration, save current terminal type and then
# set terminal to gif and specify the output going to the file uExact.gif
#

if(output_count == 0) \
set term push;           \
set terminal gif medium animate transparent opt delay 10 size 320,240;                         \
set output 'uExact.gif'

#
# For each iteration load the data file as specified by the iteration count
#

if (output_count < 10)        dataFile = sprintf("uExact00%d.dat",output_count); \
else if (output_count < 100)  dataFile = sprintf("uExact0%d.dat",output_count);  \
else                          dataFile = sprintf("uExact0%d.dat",output_count)


# Set plot parameters

unset title
unset key
set yrange [-.5:1.5]

#
# Plot the data in the data file, increment iteration counter and 
# recursively call this script (using reread) to create 
# plots of subsequent data files
#
# If the output_count exceeds output_max, then 
# then reset the termainal type. 
#

if(output_count <= output_max) \
 plot dataFile w lines;                 \
 output_count = output_count + 1; \
 reread;                                \
else                                    \
set term pop                            # restore terminal type back to default 




