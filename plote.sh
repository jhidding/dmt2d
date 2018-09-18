#!/bin/bash

cat << EOF
unset multiplot
reset

set xrange [0:512]
set yrange [0:512]

set term wxt 
#background rgbcolor 'black'
#set term pdf size 12,6
#set term pngcairo background rgbcolor 'black' size 2048,1024
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

psf=1.0
#set output 'shiny.pdf'

set size square
unset colorbox
set cbrange [0:5]
set palette model HSV functions 0.10 + gray/3, gray/2, 1-gray/2
set size square
set xrange [0:256]
set yrange [0:256]
#set palette model HSV functions 0.2 + gray/3, 0.8 - gray/2, gray

plot	'$1.a3b.E.dmt.conan' w l lc rgbcolor 'white' lw 4, \
	'$1.a3b.E.dmt.conan' u 1:2:(prgb2(\$3/10.)) w l lc rgbcolor variable lw 2, \
	'$1.a3a.E.dmt.conan' w l lc rgbcolor 'white' lw 4, \
	'$1.a3a.E.dmt.conan' u 1:2:(prgb1(\$3/10.)) w l lc rgbcolor variable lw 2, \
	'$1.caustic0.dmt.conan' w l lw 1 lc rgb '#882200', \
	'$1.caustic1.dmt.conan' w l lw 1 lc rgb '#002288'

EOF
#	'a.caustic0.dmt.conan' w l lw 2 lc rgb 'white', 'a.caustic1.dmt.conan' w l lw 2 lc rgb 'white', \

