#!/bin/bash

if $2; then
gnuplot << EOF
set contour
set cntrparam levels auto 40
set table
set output 'gp_cntr.tab'
set view map
unset surface
splot '$1.dmt.conan' i 4 matrix
unset table
EOF
fi

cat << EOF
unset multiplot
reset

set xrange [50:400]
set yrange [100:450]

#set term wxt 
#background rgbcolor 'black'
#set term pdf size 10,7
#set term pngcairo background rgbcolor 'white' size 768,768
#set border linecolor rgbcolor 'yellow'
#set tics textcolor rgbcolor 'yellow'
unset key

trunc(x)=(x < 0 ? 0 : (x > 1 ? 1 : x))
rgb(r,g,b)=int(trunc(r)*255)*65536 + int(trunc(g)*255)*256 + int(trunc(b)*255)*1
irgb(r,g,b)=int(trunc(1-r)*255)*65536 + int(trunc(1-g)*255)*256 + int(trunc(1-b)*255)*1
#prgb1(x)=rgb(x/3,x,x*4./5)
#prgb2(x)=rgb(x,x*2./3,x/2)
#prgb3(x)=rgb(x,x*4./5,x/3)
prgb1(x)=rgb(1-0.5*x,1-x**(0.3),1-x**(0.3))
prgb2(x)=rgb(1-x**(0.3),1-x**(0.3),1-0.5*x)
prgb3(x)=rgb(1-x,1-x,1-x)

psf=3.0
linef=1.0
#set output 'with-vectors.pdf'
#unset xtics
#unset ytics
#set size square
unset colorbox
set cbrange [0:5]
set palette model HSV functions 0.10 + gray/3, gray/2, 1-gray/2

#set xrange [640:740]
#set yrange [530:600]
plot 	'gp_cntr.tab' u 1:2:3 w l lw 1.0*linef palette, \
	'$1.a3b.L.dmt.conan' w l lc rgbcolor 'white' lw 5*linef, \
	'$1.a3b.L.dmt.conan' u 1:2:(prgb2(\$3/10.)) w l lc rgbcolor variable lw 3*linef, \
	'$1.a3a.L.dmt.conan' w l lc rgbcolor 'white' lw 5*linef, \
	'$1.a3a.L.dmt.conan' u 1:2:(prgb1(\$3/10.)) w l lc rgbcolor variable lw 3*linef, \
	'$1.dmt.conan' i 2 w p pt 7 ps 0.5*psf lc rgb 'white', \
	'$1.dmt.conan' i 2 u 1:2:(prgb1(\$3/5.)) w p pt 7 ps 0.3*psf lc rgb var, \
	'$1.dmt.conan' i 1 w p pt 11 ps 1.0*psf lc rgb 'white', \
	'$1.dmt.conan' i 1 u 1:2:(prgb3(\$3/5.)) w p pt 11 ps 0.6*psf lc rgb var

#	'$1.cntr2.dmt.conan' w l lc rgb 'white' lw 3, \
#	'$1.cntr2.dmt.conan' w l lc rgb '#002288' lw 2, \
#	'$1.cntr1.dmt.conan' w l lc rgb '#882200' lw 2,
#	'$1.evec.dmt.conan' u (\$1-\$3/sqrt(\$3**2 + \$4**2)):(\$2-\$4/sqrt(\$3**2 + \$4**2)):(2*\$3/sqrt(\$3**2 + \$4**2)):(2*\$4/sqrt(\$3**2 + \$4**2)) with vectors nohead lc rgb '#666666' lw linef, 

#set output
EOF
