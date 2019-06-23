#############################
# Example of Horner's Method and 
# finding roots of polynomials
# Created for the YouTube channel
# https://www.youtube.com/OscarVeliz
# @author Oscar Veliz
#############################

# Long Evaluation
function longway(a, x)
	result = 0
	for i in 1:length(a)
		result += a[i] * x^(i-1)
	end
	return result
end

# Iterative Horner Evaluation
function horner(a, x)
	result = a[end]
	for i in length(a)-1:-1:1
		result = a[i] + x * result
	end
	return result
end

# Recursive Horner Evaluation
function rHorner(a, x, i)
	if i == length(a)
		return a[i]
	else
		return a[i] + x * rHorner(a,x,i+1)
	end
end

# Synthetic Division
# Given array of coefficients and value x return quotient
function quotient(a, x)
	q = zeros(length(a))
	q[end] = a[end]
	for i in length(a)-1:-1:1
		q[i] = x * q[i+1]+a[i]
	end
	return view(q, 2:length(q))#ignore remainder
end

# Differentiate using Power Rule
function ddx(a)
	d = zeros(length(a)-1)
	for i in 1:length(d)
		d[i] = a[i+1]*i
	end
	return d
end

# Find derivative at a point using Horner
function hornerDeriv(a, x)
	q = quotient(a, x)
	return horner(q,x)
end

# Newton's Method strictly on Polynomials using Horner
function newtonHorner(a, x)
	i = 0
	px = horner(a,x)
	# Power Rule
	#q = ddx(a)
	#qx = horner(q,x)
	# Ruffini's Rule
	qx = hornerDeriv(a,x)
	while abs(px) > 0.0000001 && i < 1000
		if qx == 0 # avoid divide by zero
			return NaN
		end
		x = x - px / qx
		px = horner(a,x)
		#qx = horner(q,x) # Power Rule
		qx = hornerDeriv(a,x) # Ruffini's Rule
		i = i + 1
	end
	if i == 1000 # hit the iteration cap (divierging or stuck in loop)
		return NaN
	end
	r = round(x) # check if rounding will give you a better answer
	if(abs(horner(a,r)) < abs(horner(a,x)))
		x = r
	end
	return x
end

# Assumes that there are no imaginary roots
function findRoots(a)
	roots = zeros(length(a)-1)
	for r in 1:length(roots)
		printPoly(a)
		root = newtonHorner(a,0) # initial guess of zero
		println("r",r," = ",root)
		roots[r] = root
		a = quotient(a,root) # deflation
	end
	return roots
end

# Print in polynomial form instead of array of coefficients
function printPoly(a)
	for i in length(a):-1:1
		# remove ".0" at end
		ai = string(a[i])
		if length(ai) >= 2 && ai[end] == '0' && ai[end-1] == '.'
			ai = ai[1:end-2]
		end
		if i == 1 # no x term 
			println(ai)
		elseif i == 2 # no exponent
			print(ai,"x + ")
		else
			print(ai,"x^",(i-1)," + ")
		end
	end
end

p = [-2310, 727, 382, -72, -8, 1]
println(sort(findRoots(p)))
