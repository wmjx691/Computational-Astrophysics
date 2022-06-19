nstars = 2

au = 1.49598e13


set style line 1 lc rgb 'red' pt 7
set style line 2 lc rgb 'blue' pt 7


set xrange[0:5.5]
set yrange[-1.2:1.2]
set xtics 0.5
set xlabel "Time [Year]"
set ylabel "X [AU]"

plot 'binary_002.dat' u ($2/(365*86400)):($4/au) w l ls 1 title "365day/yr",\
'binary_0020.dat' u ($2/(365*86400)):($4/au) w l ls 2 title "x2(t) by Kepler"

set term png
set out "pro1_x2.png"
replot