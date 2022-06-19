nstars = 2

au = 1.49598e13


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

set terminal x11 background rgb 'black'
set xlabel 'ylabel' tc rgb 'white'
set ylabel 'xlabel' tc rgb 'white'
set border lc rgb 'white'
set key tc rgb 'white'

set xrange[-1:1]
set yrange[-1:1]
set size square 
set xlabel "X [AU]"
set ylabel "Y [AU]"

do for [ii=1:200] {
    plot 'binary_001.dat' every ::1::ii  u ($4/au):($5/au) w l ls 1 title "",\
         'binary_001.dat' every ::ii::ii u ($4/au):($5/au) w p ls 1 title "SUN",\
         'binary_002.dat' every ::1::ii  u ($4/au):($5/au) w l ls 2 title "",\
         'binary_002.dat' every ::ii::ii u ($4/au):($5/au) w p ls 2 title "MERCURY",\
         'binary_003.dat' every ::1::ii  u ($4/au):($5/au) w l ls 3 title "",\
         'binary_003.dat' every ::ii::ii u ($4/au):($5/au) w p ls 3 title "VENUS",\
         'binary_004.dat' every ::1::ii  u ($4/au):($5/au) w l ls 4 title "",\
         'binary_004.dat' every ::ii::ii u ($4/au):($5/au) w p ls 4 title "EARTH",\
         'binary_005.dat' every ::1::ii  u ($4/au):($5/au) w l ls 5 title "",\
         'binary_005.dat' every ::ii::ii u ($4/au):($5/au) w p ls 5 title "MOON",\
         'binary_006.dat' every ::1::ii  u ($4/au):($5/au) w l ls 6 title "",\
         'binary_006.dat' every ::ii::ii u ($4/au):($5/au) w p ls 6 title "MARS",\
         'binary_007.dat' every ::1::ii  u ($4/au):($5/au) w l ls 7 title "",\
         'binary_007.dat' every ::ii::ii u ($4/au):($5/au) w p ls 7 title "JUPITER",\
         'binary_008.dat' every ::1::ii  u ($4/au):($5/au) w l ls 8 title "",\
         'binary_008.dat' every ::ii::ii u ($4/au):($5/au) w p ls 8 title "SATURN",\
         'binary_009.dat' every ::1::ii  u ($4/au):($5/au) w l ls 9 title "",\
         'binary_009.dat' every ::ii::ii u ($4/au):($5/au) w p ls 9 title "URANUS",\
         'binary_010.dat' every ::1::ii  u ($4/au):($5/au) w l ls 10 title "",\
         'binary_010.dat' every ::ii::ii u ($4/au):($5/au) w p ls 10 title "NEPTUNE",\
         'binary_011.dat' every ::1::ii  u ($4/au):($5/au) w l ls 11 title "",\
         'binary_011.dat' every ::ii::ii u ($4/au):($5/au) w p ls 11 title "PLUTO"
}

