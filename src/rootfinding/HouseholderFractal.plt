#########################
# Householder Fractal
# Created by Oscar Veliz
# youtube.com/OscarVeliz
# Based off Halley Fractal
#########################
set terminal pngcairo transparent crop size 4778, 2757 #use 3840, 2160 without crop
set output 'HouseholderFractal.png'
max = 100 #iterations, lower if taking too long (limit around 240)
complex (x, y) = x * {1, 0} + y * {0, 1}
centerx = 0
centeri = 0
r = 2
ratio = 1.7777778 #16/9 change to 1 if using square
dots = 1000 #clarity, lower if taking too long
a = complex(.5,.5) #1/2 + i/2
p(z) = z**8 + 15*z**4 - 16
dp(z) = 8*z**7 + 60*z**3
ddp(z) = 56*z**6 + 180*z**2
dddp(z) = 336*z**5 + 360*z
house(z,pz,dpz,ddpz) = z - a * (6*pz*dpz**2-3*pz**2*ddpz) / (6*dpz**3-6*pz*dpz*ddpz+pz**2*dddp(z))
householder(z, n) = n == max || abs(p(z)) < 0.0000001 ? n : householder (house(z,p(z),dp(z),ddp(z)), n + 1)
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
splot householder(complex (x, y), 0) notitle
