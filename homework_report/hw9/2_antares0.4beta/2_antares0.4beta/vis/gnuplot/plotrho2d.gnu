set xlabel "Y"
set ylabel "Y"
set size square
#set palette rgbformulae 33,13,10

xsize = 104
ysize = 104
set xrange[0:xsize]
set yrange[0:xsize]
set cbrange [0:255]

list = system('ls rXY*.ppm')

do for [fn in list] {
    plot fn binary array=(xsize,ysize) flipy format='%uchar%uchar%uchar' using ($3+$2+$1)/3. with image
    #plot fn binary array=(xsize,ysize) flipy format='%uchar' with rgbimage
    pause 0.1
}
pause -1

