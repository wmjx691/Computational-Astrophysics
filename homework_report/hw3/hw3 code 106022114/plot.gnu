nstars = 2

au = 1.49598e13


set style line 1 lc rgb 'red' pt 7
set style line 2 lc rgb 'cyan' pt 7
set style line 3 lc rgb 'yellow' pt 7

set terminal x11 background rgb 'black'
set xlabel 'ylabel' tc rgb 'white'
set ylabel 'xlabel' tc rgb 'white'
set border lc rgb 'white'
set key tc rgb 'white'

set xrange[-15:15]
set yrange[-15:15]
set size square 
set xlabel "X [AU]"
set ylabel "Y [AU]"

do for [ii=1:200] {
    plot 'binary_1.dat' every ::1::ii  u ($4/au):($5/au) w l ls 1 title "",\
         'binary_1.dat' every ::ii::ii u ($4/au):($5/au) w p ls 1 title "mass1",\
         'binary_2.dat' every ::1::ii  u ($4/au):($5/au) w l ls 2 title "",\
         'binary_2.dat' every ::ii::ii u ($4/au):($5/au) w p ls 2 title "mass2",\
         'binary_3.dat' every ::1::ii  u ($4/au):($5/au) w l ls 3 title "",\
         'binary_3.dat' every ::ii::ii u ($4/au):($5/au) w p ls 3 title "mass3"
}

