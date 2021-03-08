nstars = 2

au = 1.49598e13


set style line 1 lc rgb 'red' pt 7
set style line 2 lc rgb 'blue' pt 7

set terminal x11 background rgb 'black'
set xlabel 'ylabel' tc rgb 'white'
set ylabel 'xlabel' tc rgb 'white'
set border lc rgb 'white'
set key tc rgb 'white'

set xrange[-5:5]
set yrange[-5:5]
set size square 
set xlabel "X [AU]"
set ylabel "Y [AU]"

do for [ii=1:200] {
    plot 'binary_001.dat' every ::1::ii  u ($4/au):($5/au) w l ls 1 title "",\
         'binary_001.dat' every ::ii::ii u ($4/au):($5/au) w p ls 1 title "m1",\
         'binary_002.dat' every ::1::ii  u ($4/au):($5/au) w l ls 2 title "",\
         'binary_002.dat' every ::ii::ii u ($4/au):($5/au) w p ls 2 title "m2"
}

