set xrange[0:1]
set yrange[0:1.1]
set xlabel "x"
set ylabel "u"

do for [ii=0:1800:10] {
    fn = sprintf('advection_%05d.d',ii)
    plot fn u 2:3 w p title "Numerical", '' u 2:4 w l title "Anayltical"
    pause 0.1
}
pause -1

