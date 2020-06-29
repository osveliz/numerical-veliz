#########################
# Graph Functions in Complex Plane
# Created by Oscar Veliz
# youtube.com/OscarVeliz
# Based off Newton Fractal program
#########################
set terminal pngcairo transparent crop size 4778, 2757 #use 3840, 2160 without crop
#set terminal pngcairo size 3840, 2160
set output 'ComplexRoots.png'
max = 1000 #for the cap normalizing
complex (x, y) = x * {1, 0} + y * {0, 1}
centerx = 0
centeri = 0
r = 2 #size of of plot (distance from center)
ratio = 1.7777778 #16/9 change to 1 if using square
#ratio = 1 #16/9 change to 1 if using square
dots = 4000 #clarity, lower if taking too long
p(z) = z**3 - 1
#p(z) = z**8 + 15*z**4 - 16
cap(z) = z >= max ? 1 : z/max
set palette rgb 6, 7, 8 #color formulation
#set palette gray
unset colorbox #remove to show color box
unset border #remove to show border
unset tics #remove to show tics
set samples dots
set isosamples dots
set pm3d map
set size 1, 1 #fill allowable canvas, to create square plot use "set size square"
#set size square
set xrange [-1 * ratio * r + centerx : ratio * r + centerx]
set yrange [-1 * r + centeri : r + centeri]
splot cap(abs(p(complex (x, y)))) notitle
