/**
 * Durand-Kerner and Aberth-Ehrlich
 * for finding all roots of polynomials
 * Example created by Oscar Veliz
 * www.youtube.com/OscarVeliz
 */
package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"time"
)

var c = make(chan int, 100) //buffered channel
var roots = []complex128{}

const ZERO = complex(0, 0)
const EPS = 10e-15
const MAX = 100000000

type ruffini struct {
	p complex128
	q complex128
}

/**
 * A - B
 * @param a one array
 * @param b another array
 * @return a - b
 */
func diff(a []complex128, b []complex128) (result []complex128) {
	result = make([]complex128, len(a))
	for i := range result {
		result[i] = a[i] - b[i]
	}
	return
}

/**
 * L2 Norm |X|
 * @param x an array
 * @return |x|
 */
func norm(x []complex128) float64 {
	sum := 0.0
	for i := range x {
		sum += math.Pow(cmplx.Abs(x[i]), 2)
	}
	return math.Sqrt(sum)
}

/**
 * Horner's Method for quickly evaluating polynomials
 * @param a array of coefficients index is x power [1,2,3] = 1 + 2x + 3x^2
 * @param x input to function
 * @return p(x)
 */
func horner(a []float64, x complex128) (result complex128) {
	result = complex(a[len(a)-1], 0)
	for i := len(a) - 2; i >= 0; i-- {
		result = complex(a[i], 0) + x*result
	}
	return
}

/**
 * Evaluate the polynomial and derivative
 * @param a array of coefficients index is x power [1,2,3] = 1 + 2x + 3x^2
 * @param x input to function
 * @return {p(x),q(x)}
 */
func hornerruffini(a []float64, x complex128) ruffini {
	n := len(a)
	q := complex(0, 0)
	p := complex(a[n-1], 0)
	for i := n - 2; i >= 0; i-- {
		q = q*x + p
		p = complex(a[i], 0) + x*p
	}
	return ruffini{p, q}
}

/**
 * Durand-Kerner
 * @param a array of coefficients index is x power [1,2,3] = 1 + 2x + 3x^2
 * @param r array of roots to iterate on
 * @param index estimated root index
 */
func durandkerner(a []float64, r []complex128, index int) {
	denominator := complex(1, 0)
	for i := 0; i < len(r); i++ {
		if i != index {
			denominator *= (r[index] - r[i])
		}
	}
	roots[index] = r[index] - horner(a, r[index])/denominator
	c <- index
}

/**
 * Aberth-Ehrlich
 * @param a array of coefficients index is x power [1,2,3] = 1 + 2x + 3x^2
 * @param r array of roots to iterate on
 * @param index estimated root index
 */
func aberthehrlich(a []float64, r []complex128, index int) {
	pofz := hornerruffini(a, r[index])
	newt := pofz.p / pofz.q
	beta := complex(0, 0)
	for i := 0; i < len(r); i++ {
		if i != index {
			beta += 1.0 / (r[index] - r[i])
		}
	}
	roots[index] = roots[index] - newt/(1.0-newt*beta)
	c <- index
}

/**
 * Durand-Kerner
 * @param a array of coefficients index is x power [1,2,3] = 1 + 2x + 3x^2
 * @param durand true when using durand-kerner, false for aberth-ehrlich
 * @return array of complex roots
 */
func solve(a []float64, durand bool) []complex128 {
	if durand {
		fmt.Println("Durand-Kerner")
	} else {
		fmt.Println("Aberth-Ehrlich")
	}
	//initialize starting circle
	roots = make([]complex128, len(a)-1)
	n := len(roots)
	img := complex(0.0, 1.0)
	rad := math.Pow(math.Abs(a[0]/a[n]), 1.0/float64(n)) //nroot(|a0/an|)
	theta := 2.0 * math.Pi / float64(n)                  //evenly space circle
	offset := theta / float64(n+1)                       //don't start with real number
	for k := 0; k < n; k++ {
		roots[k] = complex(rad, 0.0) * cmplx.Exp(img*complex(theta*float64(k)+offset, 0.0))
	}
	//begin root finding
	old := make([]complex128, n)
	loop := 0
	for loop = 0; norm(diff(old, roots)) > EPS && loop < MAX; loop++ {
		//fmt.Println(loop, roots)//uncomment to show iterations
		copy(old, roots)
		if durand {
			for i := range roots { //iterate in background
				go durandkerner(a, old, i)
			}
		} else {
			for i := range roots {
				go aberthehrlich(a, old, i)
			}
		}
		for i := 0; i < len(roots); i++ { //wait for them to finish
			<-c
		}
	}
	//fmt.Println(loop, roots)//uncomment to show last iteration
	fmt.Println(loop, "iterations")
	fmt.Println(norm(diff(old, roots)), "last step")
	return roots
}

/**
 * Run both solvers
 * @param c coefficients with index corresponding to x-power
 */
func compare(c []float64) {
	fmt.Println(c) //uncomment to show last iteration
	start := time.Now()
	fmt.Println(solve(c, true))
	fmt.Printf("%s\n", time.Since(start))
	start = time.Now()
	fmt.Println(solve(c, false))
	fmt.Printf("%s\n\n", time.Since(start))
}

/**
 * Main
 */
func main() {
	compare([]float64{-1, 0, 0, 1})                  //x^3 - 1
	compare([]float64{-16, 0, 0, 0, 15, 0, 0, 0, 1}) //x^8 + 15x^4 - 16
	compare([]float64{-2310, 727, 382, -72, -8, 1})  //(x-2)(x+3)(x-5)(x+7)(x-11)
	compare([]float64{12, -28, 23, -8, 1})           //(x-1)(x-3)(x-2)^2
	compare([]float64{16, -32, 24, -8, 1})           //(x-2)^4
	//w(x) = (x-1)(x-2)...(x-20) below
	//compare([]float64{2432902008176640000,-8752948036761600000,13803759753640704000,-12870931245150988800,8037811822645051776,-3599979517947607200,1206647803780373360,-311333643161390640,63030812099294896,-10142299865511450,1307535010540395,-135585182899530,11310276995381,-756111184500,40171771630,-1672280820,53327946,-1256850,20615,-210,1})
}
