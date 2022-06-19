#防止記憶體崩毀
set terminal x11

set style line 1 lc rgb 'blue' pt 7


#set xrange[0:1.3e6]
#set yrange[0:1e-5]
set size square 
set xlabel "R [cm]"
set ylabel "rho [/cm^3]"

plot 'nbody_0 1.dat' u ($2):($3) w l ls 1 title "rho"

set term png
set out "1.png"
replot