
au = 1.49598e13

#防止記憶體崩毀
set terminal x11

set title "Time Step for 0.1 yrs"
set xrange[-30:30]
set yrange[-30:30]
set size square 
set xlabel "X [AU]"
set ylabel "Y [AU]"

plot 'binary_0011.dat' u ($4/au):($5/au) w l lc 'red' title "mass1",\
'binary_0021.dat' u ($4/au):($5/au) w l lc 'blue' title "mass2",\
'binary_0031.dat' u ($4/au):($5/au) w l lc 'green' title "mass3"

set term png
set out "pro3_a_0.1.png"
replot