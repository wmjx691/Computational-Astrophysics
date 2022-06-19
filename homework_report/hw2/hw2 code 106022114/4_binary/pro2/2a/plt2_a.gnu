nstars = 2

au = 1.49598e13


set style line 1 lc rgb 'red' pt 7
set style line 2 lc rgb 'blue' pt 7

set xlabel 'ylabel'
set ylabel 'xlabel'

set xrange[-4:4]
set yrange[-4:4]
set size square 
set xlabel "X [AU]"
set ylabel "Y [AU]"

plot 'binary_001.dat' u ($4/au):($5/au) w p ls 1 title "mass1",\
'binary_002.dat' u ($4/au):($5/au) w p ls 2 title "mass2"

set term png
set out "pro2_a.png"
replot