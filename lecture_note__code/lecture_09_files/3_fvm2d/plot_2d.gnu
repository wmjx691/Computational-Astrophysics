set xrange[0:1]
set yrange[0:1]
set zrange[0:1]
set xlabel "X"
set ylabel "Y"
set size square

interval = 5

set pm3d map
set palette rgbformulae 33,13,10

do for [ii=0:800:interval] {
    fn = sprintf('advection_%05d.d',ii)
    splot fn u 1:2:3 
}
pause -1

