###########################
# Jarratt Fractal Creator
# Created by Oscar Veliz
# youtube.com/OscarVeliz
# Based off Newton Fractal
###########################
set terminal pngcairo transparent crop size 4778, 2757 #use 3840, 2160 without crop
set output 'JarrattFractal.png'
max = 100 #iterations, lower if taking too long
complex (x, y) = x * {1, 0} + y * {0, 1}
dist (x, y) = sqrt(x*x+y*y)
centerx = 0
centeri = 0
r = 2 #size of of plot (distance from center)
ratio = 1.7777778 #16/9 change to 1 if using square
dots = 1000 #clarity, lower if taking too long
a = complex(0.5, 0.5) #1/2+i/2
eps = 0.0000001
f(z) = z**4/4 - z
jar2(znm2,znm1,z,fz2,fz1,fz) =  z + a*( 0.5*((znm1-z)**2.0*(fz-fz2)+(znm2-z)**2.0*(fz1-fz))/((znm1-z)*(fz-fz2)+(znm2-z)*(fz1-fz)))
jar1(znm2,znm1,z,fz2,fz1,fz,n) = n == max || abs(z-znm1) < eps ? n : jar1(znm1,z,jar2(znm2,znm1,z,fz2,fz1,fz) ,fz1,fz,f(jar2(znm2,znm1,z,fz2,fz1,fz)),n+1)
jar0(znm2,znm1,z) = jar1(znm2,znm1,z,f(znm2),f(znm1),f(z),0)
jar(z,theta) = z == {0,0} ? jar0(-0.001,0.001,0) : jar0(abs(z)*1.001*exp({0,1}*theta),abs(z)*0.999*exp({0,1}*theta),z)
jarratt(x,y) = x == 0 && y > 0 ? jar(complex(x,y),pi/2) : x == 0 ? jar(complex(x,y),3*pi/2) : jar(complex(x,y),atan(y/x)+ (x<0 ? pi : 0));
set palette rgb 23, 28, 3 #color formulation
unset colorbox #remove to show color box
unset border #remove to show border
unset tics #remove to show tics
set samples dots
set isosamples dots
set pm3d map
set size 1, 1 #fill allowable canvas, to create square plot use "set size square"
set xrange [-1 * ratio * r + centerx : ratio * r + centerx]
set yrange [-1 * r + centeri : r + centeri]
splot jarratt(x, y) notitle
