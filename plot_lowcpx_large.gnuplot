# Gnuplot script file for plotting data in file "force.dat"
# This file is called   force.p
set autoscale                        # scale axes automatically
set key left top
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Different rulesets"
set xlabel "Time steps"
set ylabel "Compressed length"
plot  "data/out82.dat" using 1:2 title 'Rule 82 - gzip' with linespoints, \
      "data/out88.dat" using 1:2 title 'Rule 88 - gzip' with linespoints, \
      "data/out167.dat" using 1:2 title 'Rule 167 - gzip' with linespoints
set terminal pop
set output
replot
