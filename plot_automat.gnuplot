reset
set term pngcairo size 2000,2000
set output ARG2
unset border
unset xtics
unset ytics
unset colorbox
set yrange [ARG3:0] reverse
set autoscale
set palette grey negative
plot ARG1 matrix with image
replot
