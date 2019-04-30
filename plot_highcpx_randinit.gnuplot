# Gnuplot script file for plotting data in file "force.dat"
# This file is called   force.p
set terminal push
set terminal pngcairo
set output 'figures/rules_highcpx_800bits800ts_rand.png'
set   autoscale                        # scale axes automatically
set key left top
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Different rulesets"
set xlabel "Time steps"
set ylabel "Compressed length"
plot  "data/out101.dat" using 1:2 title 'Rule 101 - gzip' with linespoints, \
      "data/out30.dat" using 1:2 title 'Rule 30 - gzip' with linespoints, \
      "data/out89.dat" using 1:2 title 'Rule 89 - gzip' with linespoints, \
      "data/out101_paq.dat" using 1:2 title 'Rule 101 - paq' with linespoints, \
      "data/out30_paq.dat" using 1:2 title 'Rule 30 - paq' with linespoints, \
      "data/out89_paq.dat" using 1:2 title 'Rule 89 - paq' with linespoints
set terminal pdf size 6, 5
set output "figures/rules_highcpx_800bits800ts_rand.pdf"
replot
set terminal pop
set output
replot