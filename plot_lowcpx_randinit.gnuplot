# Gnuplot script file for plotting data in file "force.dat"
# This file is called   force.p
set terminal push
set terminal pngcairo
set output 'figures/rules_lowcpx_800bits800ts_rand.png'
set   autoscale                        # scale axes automatically
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
      "data/out167.dat" using 1:2 title 'Rule 167 - gzip' with linespoints, \
      "data/out82_paq.dat" using 1:2 title 'Rule 82 - paq' with linespoints, \
      "data/out88_paq.dat" using 1:2 title 'Rule 88 - paq' with linespoints, \
      "data/out167_paq.dat" using 1:2 title 'Rule 167 - paq' with linespoints
set terminal pdf size 6, 5
set output "figures/rules_lowcpx_800bits800ts_rand.pdf"
replot
set terminal pop
set output
replot