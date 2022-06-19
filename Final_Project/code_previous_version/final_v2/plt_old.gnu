
au = 1.49598e13

#防止記憶體崩毀
set terminal x11

set title "Time step for 0.01 yrs in 10 yrs"
set xrange[-10:10]
set yrange[-10:10]
set size square 
set xlabel "X [AU]"
set ylabel "Y [AU]"

plot 'binary_1.dat' u ($4/au):($5/au) w l lc 'red' title "mass1",\
'binary_2.dat' u ($4/au):($5/au) w l lc 'blue' title "mass2",\
'binary_3.dat' u ($4/au):($5/au) w l lc 'green' title "mass3"

set term png
set out "pro1_result.png"
replot