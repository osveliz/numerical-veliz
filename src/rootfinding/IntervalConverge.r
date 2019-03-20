############################
# Demonstration of the Interval of Convergence
# as part of online learning video series on
# YouTube created by Oscar Veliz.
# https://www.youtube.com/OscarVeliz
# @author Oscar Veliz
############################
f <- function(x) atan(x)
fp <- function(x) 1.0 / (x^2 + 1)
plotNewton <- function(xn = 1, gam = 1, max = 5){
	pdf(paste("newton", xn, "gamma", gam, "max", max, ".pdf", sep = ""))
	plot(f, -4, 4, xlab = paste("start =", xn,"gamma =",gam, sep = " "), ylab = " ", cex.lab = 1.4, cex.axis = 1.4,)
	title(main = "", sub = "x", cex.sub = 1.4, cex.main = 1.5 , xlab = " ", ylab = "arctan(x)", cex.axis = 1.4, cex.lab = 1.4, line = 2)
	title(main = expression(paste("Newton's Method with ", gamma)), cex.main = 1.5 , xlab = "", ylab = "", line = .7)
	abline(a = 0, b = 0)#x axis
	colors <- rainbow(max)
	colors[2] <- c('darkgoldenrod1')#otherwise you can't see it
	colors[3] <- c('chartreuse3')
	for(n in 1:max){
		segments(xn, 0 , xn, f(xn), lty = "dashed", col = colors[n])#connect to axis
		points(xn, f(xn), col = colors[n])#draw point (xn, f(xn))
		text(xn, y = 0.05, labels = bquote(x[.(n)]), col = colors[n])#label x's
		abline(a = (fp(xn) * (-xn) + f(xn)), b = fp(xn), col = colors[n])#y=a+bx
		xn = xn - gam * f(xn) / fp(xn)#newton's method
	}
}
computeInterval <- function(ld = -1.5, lc = -1, rc = 1, rd = 1.5, eps = 10^-7, gam = 1){
	pdf(paste("interval", ld, lc, rc, rd, gam, ".pdf", sep = "~"))
	while(abs(lc - ld) >= eps || abs(rd - rc) >= eps){
		print(paste(ld,"[",lc,",",rc,"]",rd))
		if(abs(lc - ld) >= eps){#left side
			converging = FALSE
			x = (lc + ld) / 2.0 #midpoint between conv and div
			for(n in 1:4){#test a few iterations to see if conv
				x = x - gam * f(x) / fp(x)
				if(x >= lc && x <= rc)
					converging = TRUE
			}
			if (converging)
				lc = (lc + ld) / 2.0
			else
				ld = (lc + ld) / 2.0
		}
		if(abs(rc - rd) >= eps){#repeat with right side
			converging = FALSE
			x = (rc + rd) / 2.0
			for(n in 1:4){
				x = x - gam * f(x) / fp(x)
				if(x >= lc && x <= rc)
					converging = TRUE
			}
			if (converging)
				rc = (rc + rd) / 2.0
			else
				rd = (rc + rd) /2.0
		}
	}
	print(paste(ld,"[",lc,",",rc,"]",rd))
	#plot the interval
	plot(f, -4, 4, xlab = paste("gamma =", gam, "[", lc, ",", rc,"]", sep = " "), ylab = " ", cex.lab = 1.4, cex.axis = 1.4,)
	title(main = "", sub = "x", cex.sub = 1.4, cex.main = 1.5 , xlab = " ", ylab = "arctan(x)", cex.axis = 1.4, cex.lab = 1.4, line = 2)
	title(main = expression(paste("Interval of Convergence with ", gamma)), cex.main = 1.5 , xlab = "", ylab = "", line = .7)
	abline(a = 0, b = 0)#x axis
	for(x in seq(-4, 4, by = 0.1)){
		if(x >= lc && x <= rc)
			points(x, 0, col = c('chartreuse3'), pch = 15)
		else
			points(x, 0, col = c('firebrick3'), pch = 4)		
	}
}
plotNewton(xn = 4, gam = .1)
computeInterval(ld = -2, rd = 2, gam = .9)
