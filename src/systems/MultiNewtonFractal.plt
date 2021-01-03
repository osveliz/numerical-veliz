#########################
# Multivariate Newton Fractal Creator
# Based off Newton Fractal program
# Created by Oscar Veliz
# youtube.com/OscarVeliz
#########################
set terminal pngcairo transparent crop size 4755, 2718 #use 3840, 2160 without crop
set output 'MultiVarNewtonFractal.png'
max = 100 #iterations, lower if taking too long (limit around 240)
complex (x, y) = x * {1, 0} + y * {0, 1}
centerx = 0
centery = 0
r = 2 #size of of plot (distance from center)
ratio = 1.7777778 #16/9 change to 1 if using square
dots = 1000 #clarity, lower if taking too long
a = 1
#x^2-y-1=0
#x-y^2+1=0
cx(x,y) = x - a * ( (2*y)*(x**2-y-1)/(4*x*y-1) - (x-y**2+1)/(4*x*y-1) )
cy(x,y) = y - a * ( (-2*x*(x-y**2+1))/(4*x*y-1) + (x**2-y-1)/(4*x*y-1) )
l2norm(x1, y1, x2, y2) = sqrt( (x2-x1)**2 + (y2-y1)**2 )
eps = 0.000001
newt(xn, yn, xnm1, ynm1, n) = n == max || l2norm(xn,yn,xnm1,ynm1) < eps ? n : newt(cx(xn,yn),cy(xn,yn),xn,yn, n+1)
newton(x0,y0) = newt(cx(x0,y0) ,cy(x0,y0) ,x0 ,y0 ,0)
set palette rgb 23, 28, 3 #color formulation
unset colorbox #remove to show color box
unset border #remove to show border
unset tics #remove to show tics
set samples dots
set isosamples dots
set pm3d map
set size 1, 1 #fill allowable canvas, to create square plot use "set size square"
set xrange [-1 * ratio * r + centerx : ratio * r + centerx]
set yrange [-1 * r + centery : r + centery]
splot newton(x, y) notitle