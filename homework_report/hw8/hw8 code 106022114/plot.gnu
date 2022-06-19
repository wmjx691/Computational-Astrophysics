set terminal x11

set xrange[0:1]
set yrange[0:1.1]
set xlabel "x"
set ylabel "u"

set style line 1 lc rgb 'red' pt 7 lw 2.5
set style line 2 lc rgb 'orange' pt 7 lw 2.5
set style line 3 lc rgb 'yellow' pt 7 lw 2.5
set style line 4 lc rgb 'green' pt 7 lw 2.5
set style line 5 lc rgb 'blue' pt 7 lw 2.5
set style line 6 lc rgb 'violet' pt 7 lw 2.5


plot 'advection_00000.d' u 2:3 w l ls 1 title "t=0.000(s)",\
     'advection_00010.d' u 2:3 w l ls 2 title "t=0.013(s)",\
     'advection_00020.d' u 2:3 w l ls 3 title "t=0.026(s)",\
     'advection_00030.d' u 2:3 w l ls 4 title "t=0.039(s)",\
     'advection_00040.d' u 2:3 w l ls 5 title "t=0.052(s)",\
     'advection_00046.d' u 2:3 w l ls 6 title "t=0.060(s)"
  
set term png
set out "t13_4.png"
replot

