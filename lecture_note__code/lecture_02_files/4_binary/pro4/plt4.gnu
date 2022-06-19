au = 1.49598e13

#防止記憶體崩毀
set terminal x11

set style line 1 lc rgb 'red' pt 7
set style line 2 lc rgb 'cyan' pt 7
set style line 3 lc rgb 'yellow' pt 7
set style line 4 lc rgb 'blue' pt 7
set style line 5 lc rgb 'gray' pt 7
set style line 6 lc rgb 'orange' pt 7
set style line 7 lc rgb 'brown' pt 7
set style line 8 lc rgb 'khaki' pt 7
set style line 9 lc rgb 'cyan' pt 7
set style line 10 lc rgb 'web-blue' pt 7
set style line 11 lc rgb 'gray' pt 7

set xrange[-5:5]
set yrange[-5:5]
set size square 
set xlabel "X [AU]"
set ylabel "Y [AU]"

plot 'binary_001.dat' u ($4/au):($5/au) w l ls 1 title "SUN",\
        'binary_002.dat' u ($4/au):($5/au) w l ls 2 title "MERCURY",\
        'binary_003.dat' u ($4/au):($5/au) w l ls 3 title "VENUS",\
        'binary_004.dat' u ($4/au):($5/au) w l ls 4 title "EARTH",\
        'binary_005.dat' u ($4/au):($5/au) w l ls 5 title "MOON",\
        'binary_006.dat' u ($4/au):($5/au) w l ls 6 title "MARS",\
        'binary_007.dat' u ($4/au):($5/au) w l ls 7 title "JUPITER",\
        'binary_008.dat' u ($4/au):($5/au) w l ls 8 title "SATURN",\
        'binary_009.dat' u ($4/au):($5/au) w l ls 9 title "URANUS",\
        'binary_010.dat' u ($4/au):($5/au) w l ls 10 title "NEPTUNE",\
        'binary_011.dat' u ($4/au):($5/au) w l ls 11 title "PLUTO",\

set term png
set out "pro42.png"
replot
