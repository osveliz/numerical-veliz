############################
# Visualization of Newton's and Halley's Method
# demo of the Convergance Interval for both
# https://www.youtube.com/OscarVeliz
# @author Oscar Veliz
############################
hyperx <- a <- b <- c <- 0
center <- 0
range <- 4
cap <- 2000 #maximum value for hyperbola
stepsize <- .001 #when generating points
eps <- 10^-7 #epsilon
hyperbola <- function(x) (x - hyperx + c) / (a*(x - hyperx)+b)
max <- 5
colors <- rainbow(max)
colors[2] <- c('darkgoldenrod1')#otherwise you can't see it
colors[3] <- c('chartreuse3')#same here	
############################
# Change the following functions
############################
name <- "arctan(x)"
f <- function(x) atan(x)
fp <- function(x) 1.0 / (x^2 + 1)
fp2 <- function(x) (-2*x) / ((x^2+1)^2)
newtstep <- function(x) f(x) / fp(x)
halleystep <- function(x){
	fx = f(x)
	fpx = fp(x)
	(2*fx*fpx) / (2*fpx^2 - fx*fp2(x))
}
#####################
# Plots Newton's or Halley's Method
# @param xn starting x
# @param gam gamma (damping value)
# @param useNewt TRUE Newton, FALSE Halley
##################################
plotNewton <- function(xn = 1, gam = 1, useNewt = TRUE){
	if(useNewt)
		method = "Newton"
	else
		method = "Halley"
	print(paste(method,"Method"))
	svg(paste(method,"Plot.svg",sep = ""), bg = 'transparent')
	plot(f, center - range, center + range, xlab = paste("start =", xn,"gamma =",gam, sep = " "), ylab = " ", cex.lab = 1.4, cex.axis = 1.4, bg = 'transparent')
	title(main = "", sub = "x", cex.sub = 1.4, cex.main = 1.5 , xlab = " ", ylab = name, cex.axis = 1.4, cex.lab = 1.4, line = 2)
	title(main = paste(method,"'s Method", sep = ""), cex.main = 1.5 , xlab = "", ylab = "", line = .7)
	abline(a = 0, b = 0)#x axis
	if(useNewt){#Newton Plot
		for(n in 1:max){
			print(paste("x",n," = ",xn,sep=""))
			segments(xn, 0 , xn, f(xn), lty = "dashed", col = colors[n])#connect to axis
			points(xn, f(xn), col = colors[n])#draw point (xn, f(xn))
			text(xn, y = 0.05, labels = bquote(x[.(n)]), col = colors[n])#label x's
			abline(a = (fp(xn) * (-xn) + f(xn)), b = fp(xn), col = colors[n])#y=a+bx
			xn = xn - gam * newtstep(xn)#newton's method
		}
	}
	else{#Halley Plot
		xpoints = seq(center - range, center + range, stepsize)
		for(n in 1:max){
			print(paste("x",n," = ",xn,sep=""))
			denom = 2*fp(xn)^2-f(xn)*fp2(xn)
			hyperx <<- xn
			a <<- (-fp2(xn)) / denom
			b <<- (2*fp(xn)) / denom
			c <<- (2*f(xn)*fp(xn)) / denom
			segments(xn, 0 , xn, f(xn), lty = "dashed", col = colors[n])#connect to axis
			points(xn, f(xn), col = colors[n], pch = '.' )#draw point (xn, f(xn))
			text(xn, y = 0.05, labels = bquote(x[.(n)]), col = colors[n])#label x's
			ypoints = hyperbola(xpoints)
			#sanitize to remove asymptote connection
			for(i in 1:length(ypoints))
				if(abs(ypoints[i])>cap)
					ypoints[i] = NA
			lines(xpoints,ypoints,col = colors[n])
			points(xn, f(xn), col = colors[n])#draw point (xn, f(xn))
			xn = xn - gam * halleystep(xn)#halley's method
			if(abs(f(xn))<eps)
				break;
		}
	}
	print(paste("x = ",xn,sep=""))
}
#####################
# Plots Halley's Method
# @param xn starting x
# @param gam gamma (damping value)
# @param max maximum number of iterations to display
##################################
plotHalley <- function(xn = 1, gam = 1) plotNewton(xn,gam,useNewt=FALSE)

#####################
# Plots Convergance Interval
# Default values assume converges at center and
# diverges at ends. Otherwise makes incorrect plot.
# @param ld left diverge value
# @param lc left converge value
# @param rc right converge value
# @param rd right diverge value
# @param gam gamma damping value
# @param useNewt TRUE Newton, FALSE Halley
##################################
plotInterval <- function(ld = center - range, lc = center, rc = center, rd = center + range, gam = 1, useNewt = TRUE){
	if(useNewt)
		method = "Newton"
	else
		method = "Halley"
	print(paste(method,"Convergence Interval"))
	svg(paste(method,"ConvergenceInterval",".svg", sep = ""), bg = 'transparent')
	while(abs(lc - ld) >= eps || abs(rd - rc) >= eps){
		print(paste(ld,"[",lc,",",rc,"]",rd))
		if(abs(lc - ld) >= eps){#left side
			converging = FALSE
			x = (lc + ld) / 2.0 #midpoint between conv and div
			for(n in 1:10){#test a few iterations to see if conv
				if(useNewt)
					x = x - gam * newtstep(x)
				else
					x = x - gam * halleystep(x)
				if(x >= lc && x <= rc)
					converging = TRUE
				if(x < ld || x > rd)
					break
			}
			if (converging)
				lc = (lc + ld) / 2.0
			else
				ld = (lc + ld) / 2.0
		}
		if(abs(rc - rd) >= eps){#repeat with right side
			converging = FALSE
			x = (rc + rd) / 2.0
			for(n in 1:10){
				if(useNewt)
					x = x - gam * newtstep(x)
				else
					x = x - gam * halleystep(x)
				if(x >= lc && x <= rc)
					converging = TRUE
				if(x < ld || x > rd)
					break
			}
			if (converging)
				rc = (rc + rd) / 2.0
			else
				rd = (rc + rd) /2.0
		}
	}
	print(paste(ld,"[",lc,",",rc,"]",rd))
	#plot the interval
	plot(f, center - range, center + range, xlab = paste("gamma =", gam, "[", lc, ",", rc,"]", sep = " "), ylab = " ", cex.lab = 1.4, cex.axis = 1.4,)
	title(main = "", sub = "x", cex.sub = 1.4, cex.main = 1.5 , xlab = " ", ylab = name, cex.axis = 1.4, cex.lab = 1.4, line = 2)
	title(main = paste(method, "Convergence Interval",sep = " "), cex.main = 1.5 , xlab = "", ylab = "", line = .7)
	abline(a = 0, b = 0)#x axis
	for(x in seq(center - range + 0.1, center + range - 0.1, by = range/40))
		if(x >= lc && x <= rc)
			points(x, 0, col = c('chartreuse3'), pch = 15)
		else
			points(x, 0, col = c('firebrick3'), pch = 4)
}

### main
start <- 1.5
plotNewton(start)
plotHalley(start)
plotInterval()
plotInterval(useNewt = FALSE)
