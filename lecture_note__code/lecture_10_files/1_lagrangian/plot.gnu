#set xrange[0:1]
set ylabel "rho, u, P"

list = system('ls hydro*.d')

#a = 2 # m
#set xlabel "M"

a = 3 # r/rmax
set xlabel "r/rmax"

do for [fn in list] {
    #fn = sprintf('hydro_%05d.d',ii)
    plot fn u a:5 w lp title "rho", '' u a:4 w lp title "v", '' u a:($6/10) w lp title "P/10"
    pause 0.1
}
pause -1

