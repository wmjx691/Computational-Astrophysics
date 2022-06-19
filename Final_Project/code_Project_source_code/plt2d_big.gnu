au = 1.49598e8


#set style line 1 lc rgb 'gray' pt 7
#set style line 2 lc rgb 'gray' pt 7
set style line 3 lc rgb 'green' pt 7
#set style line 4 lc rgb 'gray' pt 7
set style line 5 lc rgb 'orange' pt 7
#set style line 6 lc rgb 'brown' pt 7
#set style line 7 lc rgb 'green' pt 7
#set style line 8 lc rgb 'cyan' pt 7
#set style line 9 lc rgb 'gray' pt 7
set style line 10 lc rgb 'red' pt 7
set style line 11 lc rgb 'white' pt 7

set term gif animate background rgb 'black'
set term gif animate delay 40
set output 'output2d.gif'

set datafile separator ","
#set terminal x11 background rgb 'black'
set xlabel 'xlabel' tc rgb 'white'
set ylabel 'ylabel' tc rgb 'white'
#set zlabel 'zlabel' tc rgb 'white'
set border lc rgb 'white'
set key tc rgb 'white'

set title tc rgb 'white'
set title 'Check the trajectory from Earth to Jupiter 2D'


set xrange[-4.5:-4.45]
set yrange[2.6:2.7]
#set zrange[-0.5:0.5]
set size square 
set xlabel "X [AU]"
set ylabel "Y [AU]"
#set zlabel "Z [AU]"


do for [ii=785:800:1] {

    plot 'planet_10.dat' every ::1::ii  u ($5/au):($6/au):($7/au) w l ls 10 title "",\
         'planet_10.dat' every ::ii::ii u ($5/au):($6/au):($7/au) w p ls 10 title "Sun",\
         'planet_3.dat' every ::1::ii  u ($5/au):($6/au):($7/au) w l ls 3 title "",\
         'planet_3.dat' every ::ii::ii u ($5/au):($6/au):($7/au) w p ls 3 title "Earth",\
         'planet_5.dat' every ::1::ii  u ($5/au):($6/au):($7/au) w l ls 5 title "",\
         'planet_5.dat' every ::ii::ii u ($5/au):($6/au):($7/au) w p ls 5 title "Jupiter",\
         'stell.dat' every ::1::ii  u ($5/au):($6/au):($7/au) w l ls 11 title "",\
         'stell.dat' every ::ii::ii u ($5/au):($6/au):($7/au) w p ls 11 title "stellite"
         
}


#         'planet_1.dat' every ::1::ii  u ($5/au):($6/au):($7/au) w l ls 1 title "",\
#         'planet_1.dat' every ::ii::ii u ($5/au):($6/au):($7/au) w p ls 1 title "Mercury",\
#         'planet_2.dat' every ::1::ii  u ($5/au):($6/au):($7/au) w l ls 2 title "",\
#         'planet_2.dat' every ::ii::ii u ($5/au):($6/au):($7/au) w p ls 2 title "Venus",
#         'planet_4.dat' every ::1::ii  u ($5/au):($6/au):($7/au) w l ls 4 title "",\
#         'planet_4.dat' every ::ii::ii u ($5/au):($6/au):($7/au) w p ls 4 title "Mars",
#         'planet_6.dat' every ::1::ii  u ($5/au):($6/au):($7/au) w l ls 6 title "",\
#         'planet_6.dat' every ::ii::ii u ($5/au):($6/au):($7/au) w p ls 6 title "Saturn",\
#         'planet_7.dat' every ::1::ii  u ($5/au):($6/au):($7/au) w l ls 7 title "",\
#         'planet_7.dat' every ::ii::ii u ($5/au):($6/au):($7/au) w p ls 7 title "Uranus",\
#         'planet_8.dat' every ::1::ii  u ($5/au):($6/au):($7/au) w l ls 8 title "",\
#         'planet_8.dat' every ::ii::ii u ($5/au):($6/au):($7/au) w p ls 8 title "Neptune",
#         'planet_9.dat' every ::1::ii  u ($5/au):($6/au):($7/au) w l ls 9 title "",\
#         'planet_9.dat' every ::ii::ii u ($5/au):($6/au):($7/au) w p ls 9 title "Pluto"