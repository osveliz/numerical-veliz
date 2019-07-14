###########################
# Laguerre Fractal Creator
# Created by Oscar Veliz
# youtube.com/OscarVeliz
# Based off Newton Fractal
###########################
set terminal pngcairo transparent crop size 4778, 2757 #use 3840, 2160 without crop
set output 'LaguerreFractal.png'
max = 100 #iterations, lower if taking too long (limit around 240)
centerx = 0
centeri = 0
r = 2 #size of of plot (distance from center)
ratio = 1.7777778 #16/9 change to 1 if using square
dots = 1000 #clarity, lower if taking too long
a = {.5,.5} #1/2+i/2
complex (x, y) = x * {1, 0} + y * {0, 1}
p(z) = z**8 + 15*z**4 - 16
dp(z) = 8*z**7 + 60*z**3
ddp(z) = 56*z**6 + 180*z**2
deg = 8
maxabs(u,v) = abs(u) > abs(v) ? u : v
lag3(z,G,sr) = z - a * deg / maxabs(G + sr,G - sr)
lag2(z,G,H) = lag3(z,G,sqrt((deg-1)*(deg*H-G**2)))
lag1(z,G) = lag2(z,G,G**2 - ddp(z)/p(z))
laguerre(z, n) = n == max || abs(p(z)) < 0.00000001 ? n : dp(z)==0 && ddp(z)==0 ? max : laguerre (lag1(z,dp(z)/p(z)), n + 1)
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
splot laguerre(complex (x, y), 0) notitle
