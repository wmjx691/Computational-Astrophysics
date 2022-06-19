nstars = 2

au = 1.49598e13
a = 0.5*(au+0.5349192658994514E+14)
e = au/a-1

G= 6.67428e-8 
msun = 1.989e33


set style line 1 lc rgb 'red' pt 7
set style line 2 lc rgb 'blue' pt 7

set xrange[0:5.5]
set yrange[0:1.7e13]
set xlabel "Time [Year]"
set ylabel "(vy2)^2 [cm^2/s^2]"

plot 'binary_002.dat' u ($2/(365*86400)):(($7)**2) w l ls 2 title "mass2",\
'binary_002.dat' u ($2/(365*86400)):(((G*msun)/(a))*((1-e)/(1+e))) w l lc rgb "red" title "Velcoity at Perihelion",\
'binary_002.dat' u ($2/(365*86400)):(((G*msun)/(a))*((1+e)/(1-e))) w l lc rgb "green" title "Velcoity at Aphelion"

set term png
set out "pro2_b.png"
replot