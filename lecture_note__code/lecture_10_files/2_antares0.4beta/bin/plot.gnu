#set xrange[0:1]
#set yrange[0:1.1]
set xlabel "x"
set ylabel "u"

list = system('ls data*.tab')

#do for [ii=0:1800:10] {
do for [fn in list] {
    #fn = sprintf('den%04d.tab',ii)
    plot fn u 2:3 w lp
    pause 0.1
}
pause -1

