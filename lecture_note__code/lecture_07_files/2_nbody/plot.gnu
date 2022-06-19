nstars   = 1
noutputs = 1000
au   = 1.49598e13
rmax = 2  

set style line 1 lc rgb 'red' pt 7
set style line 2 lc rgb 'blue' pt 7
set style line 3 lc rgb 'green' pt 7

set terminal x11

#set terminal x11 background rgb 'black'
#set xlabel 'ylabel' tc rgb 'white'
#set ylabel 'xlabel' tc rgb 'white'
#set border lc rgb 'white'
#set key tc rgb 'white'

#set xrange[0:1.414213*8.0*3.1415926]
#set yrange[0.0000000:1.0000000]
set size square 
set xlabel "Time "
set ylabel "r(t) "

#do for [ii=1:noutputs] {
#    plot 'nbody_001.dat' every ::1::ii  u ($3):($4) w l ls 1 title "",\
#         'nbody_001.dat' every ::ii::ii u ($3):($4) w p ls 1 title "",\
#         'nbody_002.dat' every ::1::ii  u ($4/au):($5/au) w l ls 2 title "",\
#         'nbody_002.dat' every ::ii::ii u ($4/au):($5/au) w p ls 2 title "Mercyry",\
#         'nbody_003.dat' every ::1::ii  u ($4/au):($5/au) w l ls 3 title "",\
#         'nbody_003.dat' every ::ii::ii u ($4/au):($5/au) w p ls 3 title "Venus",\
#         'nbody_004.dat' every ::1::ii  u ($4/au):($5/au) w l ls 4 title "",\
#         'nbody_004.dat' every ::ii::ii u ($4/au):($5/au) w p ls 4 title "Earth",\
#         'nbody_005.dat' every ::1::ii  u ($4/au):($5/au) w l ls 5 title "",\
#         'nbody_005.dat' every ::ii::ii u ($4/au):($5/au) w p ls 5 title "Mars",\
#         'nbody_006.dat' every ::1::ii  u ($4/au):($5/au) w l ls 6 title "",\
#         'nbody_006.dat' every ::ii::ii u ($4/au):($5/au) w p ls 6 title "Jupyter",\
#         'nbody_007.dat' every ::1::ii  u ($4/au):($5/au) w l ls 7 title "",\
#         'nbody_007.dat' every ::ii::ii u ($4/au):($5/au) w p ls 7 title "Saturn",\
#         'nbody_008.dat' every ::1::ii  u ($4/au):($5/au) w l ls 8 title "",\
#         'nbody_008.dat' every ::ii::ii u ($4/au):($5/au) w p ls 8 title "Uranus",\
#         'nbody_009.dat' every ::1::ii  u ($4/au):($5/au) w l ls 9 title "",\
#         'nbody_009.dat' every ::ii::ii u ($4/au):($5/au) w p ls 9 title "Neptune"
#}

plot 'nbody_001.dat' u ($2):($5) w l ls 1 lc rgb "blue" title "r(t)",\
'nbody_001.dat' u ($2):(1.0) w l ls 1 lc rgb "red" title "value r=1"


set term png
set out "pro1f2.png"
replot

pause -1 



