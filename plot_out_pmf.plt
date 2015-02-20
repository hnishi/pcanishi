#set term png
#set output "out_pmf.dat.png"

set xlabel "c1"
set ylabel "c2"

set cbrange[0:3 ]
set zrange[0:3 ]

set pm3d 
#set contour
set view 0,0
#set dgrid3d 30,30
#set hidden3d
unset surface
set palette rgbformulae 33,13,10
set cntrparam levels 20
set style line 1 lt 1 #lc rgb '#000000'
#set style increment user
splot "out_pmf.dat" u 1:2:3 with pm3d

pause -1
