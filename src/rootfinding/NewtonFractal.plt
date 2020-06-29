#########################
# Newton Fractal Creator
# Created by Oscar Veliz
# youtube.com/OscarVeliz
# Based off Mandelbrot Set plotter from
# https://rosettacode.org/wiki/Mandelbrot_set#gnuplot
#########################
set terminal pngcairo transparent crop size 4778, 2757 #use 3840, 2160 without crop
set output 'NewtonFractal.png'
max = 100 #iterations, lower if taking too long (limit around 240)
complex (x, y) = x * {1, 0} + y * {0, 1}
centerx = 0
centeri = 0
r = 2 #size of of plot (distance from center)
ratio = 1.7777778 #16/9 change to 1 if using square
dots = 1000 #clarity, lower if taking too long
a = complex(0.5, 0.25) #generalization, use complex(1, 0) for normal
p(z) = z ** 3 - 1
dp(z) = 3 * z ** 2
newt(z) = z - a * p(z) / dp(z)
newton(z, n) = n == max || abs(p(z)) < 0.0000001 || dp(z) == 0 ? n : newton (newt(z), n + 1) #dp(z) == 0 condition is optional
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
splot newton(complex (x, y), 0) notitle
