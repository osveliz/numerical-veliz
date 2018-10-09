/**
 * Demonstration of Wegstein's Method
 * as part of online learning video series on
 * YouTube created by Oscar Veliz.
 * https://www.youtube.com/OscarVeliz
 * @author Oscar Veliz
 */
import kotlin.math.abs
import kotlin.math.sqrt
fun xSquaredMinusXMinusOne(x: Double): Double = x*x - x - 1
fun oneOverXMinusOne(x: Double): Double = 1.0/(x - 1.0)
fun onePlusOneOverX(x: Double): Double = 1.0 + (1.0/x)
fun sqrtXPlusOne(x: Double): Double = sqrt(x + 1.0)
fun xSquaredMinusOne(x: Double): Double = x*x -1

/**
 * Fixed Point Iteration (fpi)
 * @param x initial x
 * @param f f(x) function
 * @param g g(x) function
 * @param eps epsilon stopping criteria default 10^-7
 * @param max maximum number of iterations default 100
 * @return the value of the last completed iteration
 */
fun fpi(x: Double, f:(x: Double) -> Double, g:(x: Double) -> Double, eps: Double = 0.000001, max: Int = 100): Double{
    var i = 0
    var xn = x
    while(abs(f(xn)) > eps && i < max){
        println("$i\t$xn")
        xn = g(xn)
        i++
    }
    return xn
}

/**
 * Wegstein's method
 * @param x initial x
 * @param f f(x) function
 * @param g g(x) function
 * @param initQ initial q value default of zero will compute constant q
 * @param eps epsilon stopping criteria default 10^-7
 * @param max maximum number of iterations default 100
 * @return the value of the last completed iteration
 */
fun wegstein(x: Double, f:(x: Double) -> Double, g:(x: Double) -> Double, initQ: Double = 0.0, eps: Double = 0.000001, max: Int = 100): Double{
    var i = 0
    var xn = x
    var q = initQ
    if(initQ == 0.0){//q not specified, compute near optimal
        var x1 = x
        var x2 = g(x1)
        var x3 = g(x2)
        var a = (x3-x2) / (x2-x1)
        println("a = $a")
        q = a / (a-1)
    }
    println("q = $q")
    while(abs(f(xn)) > eps && i < max){
        println("$i\t$xn")
        xn = q*xn + (1-q)*g(xn)
        i++
    }
    return xn
}

/**
 * Wegstein's method with updating q
 * @param x initial x
 * @param f f(x) function
 * @param g g(x) function
 * @param eps epsilon stopping criteria default 10^-7
 * @param max maximum number of iterations default 100
 * @return the value of the last completed iteration
 */
fun wegstein2(start: Double, f:(x: Double) -> Double, g:(x: Double) -> Double, eps: Double = 0.000001): Double{
    var i = 1
    val xn = DoubleArray(100)
    val gn = DoubleArray(100)
    xn[0] = start
    gn[0] = g(xn[0])
    xn[1] = gn[0]
    gn[1] = g(xn[1])
    println("0\t${xn[0].toString()}")
    println("1\t${xn[1].toString()}")
    while(abs(f(xn[i])) > eps){
        xn[i+1] = (xn[i-1]*gn[i] - xn[i]*gn[i-1]) / (xn[i-1] + gn[i] - xn[i] - gn[i-1])
        i++
        println("$i\t${xn[i]}")
        gn[i] = g(xn[(i)])
    }
    return xn[i]
}

fun main(args: Array<String>) {
    // 1+1/x
    println(fpi(2.0, ::xSquaredMinusXMinusOne, ::onePlusOneOverX))
    println(wegstein(2.0, ::xSquaredMinusXMinusOne, ::onePlusOneOverX))
    println(wegstein2(2.0, ::xSquaredMinusXMinusOne, ::onePlusOneOverX))

    // 1/(x-1)
    //println(fpi(1.5, ::xSquaredMinusXMinusOne, ::oneOverXMinusOne))
    //println(wegstein(1.5, ::xSquaredMinusXMinusOne, ::oneOverXMinusOne))
    //println(wegstein2(1.5, ::xSquaredMinusXMinusOne, ::oneOverXMinusOne))

    // sqrt(x+1)
    //println(fpi(2.0, ::xSquaredMinusXMinusOne, ::sqrtXPlusOne))
    //println(wegstein(2.0, ::xSquaredMinusXMinusOne, ::sqrtXPlusOne))
    //println(wegstein2(2.0, ::xSquaredMinusXMinusOne, ::sqrtXPlusOne))

    // x^2-1
    //println(fpi(2.0, ::xSquaredMinusXMinusOne, ::xSquaredMinusOne))
    //println(wegstein(2.0, ::xSquaredMinusXMinusOne, ::xSquaredMinusOne))
    //println(wegstein2(2.0, ::xSquaredMinusXMinusOne, ::xSquaredMinusOne))
}
