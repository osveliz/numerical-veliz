{-
 - Demonstration of Dekker-Brent Method
 - as well as Bisection Method, Secant
 - Method, and Inverse Quadratic Interpolation
 - as part of online learning video series on
 - YouTube created by Oscar Veliz.
 - https://www.youtube.com/OscarVeliz
 - @author Oscar Veliz
 -}
import Data.Ratio
data Result = Result{root :: Rational , steps :: Integer}deriving (Show)  

--function f used for all exampls
f x = x^3 - x^2 - x - 1

--epsilon
eps :: Rational
eps = toRational(10**(-7))

--Biscetion Method initializer
bisect :: Rational -> Rational -> Result
bisect a b
    | a > b = bis b a 0
    | otherwise = bis a b 0
    
--Bisection body
bis :: Rational -> Rational -> Integer -> Result
bis a b i
    | abs(b - a) < eps = Result{root = c, steps = i}
    | signum(f(a)) == signum(f(c)) = bis c b (i+1)
    | otherwise = bis a c (i+1)
    where c = (a+b) / 2

--Secant Method initializer
secant :: Rational -> Rational -> Result
secant xn xm = sec xn xm 0

--Secant Method body
sec :: Rational -> Rational -> Integer -> Result
sec xn xm i
    | abs(f(xn)) < eps = Result{root = xn, steps = i}
    | otherwise = sec xnp xn (i+1)
    where xnp = xn - (f(xn)/((f(xn)-f(xm))/(xn-xm)))

--Test if x is between a and b even
--when a and b are in reverse order
between :: Rational -> Rational -> Rational -> Bool
between a b x
    | a <= b = x >= a && x <= b
    | otherwise = x >= b && x <= a

--Dekker's Method initializer
--a and b are bracketing interval
dekker:: Rational -> Rational -> Result
dekker a b
    | abs(f(a)) < abs(f(b)) = dek b a b 0
    | otherwise = dek a b a 0

--Main body of Dekker's Method
--a is interval point, opposite sign of b
--b is other interval point and current best guess
--c is previous b (on first call initially a)
--i is the iteration count (initially zero)
dek:: Rational -> Rational -> Rational -> Integer -> Result
dek a b c i
    | abs(b - a) < eps = Result{root = b, steps = i}
    | abs(f(a)) < abs(f(b)) = dek b a a i --swap not counting as iteration
    | between b m s && signum(f(s)) /= signum(f(a)) = dek a s b (i+1) --interp
    | between b m s && signum(f(s)) == signum(f(a)) = dek b s b (i+1) --interp
    | signum(f(m)) /= signum(f(a)) = dek a m b (i+1) --bisect
    | otherwise = dek b m b (i+1) --bisect
    where s = if f(b) /= f(c) then b - f(b)/((f(b)-f(c))/(b-c)) else m
          m = (a+b)/2

--Inverse Quadratic
--Compute the root of the "sideways" Lagrange polynomial given 3 points
iq :: Double -> Double -> Double -> Double
iq xn xm xm2 = xm2*(fm*fn)/((fm2-fm)*(fm2-fn))+
              xm*(fm2*fn)/((fm-fm2)*(fm-fn))+
              xn*(fm2*fm)/((fn-fm2)*(fn-fm))
              where fn = f(xn)
                    fm = f(xm)
                    fm2 = f(xm2)

--IQ Wrapper for rational numbers
iqr :: Rational -> Rational -> Rational -> Rational
iqr xn xm xm2 = toRational(iq (realToFrac(xn)) (realToFrac(xm)) (realToFrac(xm2)))

--Inverse Quadratic Interpolation using doubles
iqi :: Double -> Double -> Double -> Double
iqi xn xm xm2
    | abs(f(xn)) < realToFrac(eps) = xn
    | otherwise = iqi (iq xn xm xm2) xn xm

--Same as iqi except returning a Result
iqires :: Double -> Double -> Double -> Integer -> Result
iqires xn xm xm2 i
    | abs(f(xn)) < realToFrac(eps) = Result{root = toRational(xn), steps = i}
    | otherwise = iqires xn1 xn xm (i+1)
    where xn1 = iq xn xm xm2

--IQI Initializer using 3 points
invquadinterp xn xm xm2 = iqires xn xm xm2 0

--Brent-Dekker Method Initializer
--a and b are bracketting interval
brentdekker :: Rational -> Rational -> Result
brentdekker a b
    | abs(f(a)) < abs(f(b)) = brent b a b 0
    | otherwise = brent a b a 0

--Brent-Dekker Method with check Initializer
--Ensures Secant Method isn't causing a slowdown
--a and b are bracketting interval
brentdekkerill :: Rational -> Rational -> Result
brentdekkerill a b
    | abs(f(a)) < abs(f(b)) = brent2 b a b b False 0
    | otherwise = brent2 a b a a False 0

--Brent-Dekker main body
--a is interval point, opposite sign of b
--b is other interval point and current best guess
--c is previous b (on first call initially a)
--i is the iteration count (initially zero)
brent :: Rational -> Rational -> Rational -> Integer -> Result
brent a b c i
    | abs(b - a) <= eps || f(b) == 0 = Result{root = b, steps = i}
    | abs(f(a)) < abs(f(b)) = brent b a c i --swap not counting as iteration
    | between k b s && signum(f(s)) /= signum(f(a)) = brent a s b (i+1) --interp
    | between k b s && signum(f(s)) == signum(f(a)) = brent b s b (i+1) --interp
    | signum(f(m)) /= signum(f(a)) = brent a m b (i+1) --bisect
    | otherwise = brent b m b (i+1) --bisect
    where s = if (f(a) /= f(c) && f(b) /= f(c)) then iqr a b c else if (b-c) /= 0 || (f(b)-f(c)) /= 0 then b - f(b)/((f(b)-f(c))/(b-c)) else m
          m = (a+b)/2
          k = (3*a+b)/4

--Brent-Dekker main body with check
--a is interval point, opposite sign of b
--b is other interval point and current best guess
--c is previous b (on first call initially a)
--d is previous value of c (initially a)
--bflag true if previous iteration was bisection
--i is the iteration count (initially zero)
brent2 :: Rational -> Rational -> Rational -> Rational -> Bool-> Integer -> Result
brent2 a b c d bflag i
    | abs(b - a) <= eps || f(b) == 0 = Result{root = b, steps = i}
    | abs(f(a)) < abs(f(b)) = brent2 b a c d bflag i --swap not counting as iteration
    | test && between k b s && signum(f(s)) /= signum(f(a)) = brent2 a s b c False (i+1) --interp
    | test && between k b s && signum(f(s)) == signum(f(a)) = brent2 b s b c False (i+1) --interp
    | signum(f(m)) /= signum(f(a)) = brent2 a m b c True (i+1) --bisect
    | otherwise = brent2 b m b c True (i+1) --bisect
    where s = if (f(a) /= f(c) && f(b) /= f(c)) then iqr a b c else if (b-c) /= 0 || (f(b)-f(c)) /= 0 then b - f(b)/((f(b)-f(c))/(b-c)) else m
          m = (a+b)/2
          k = (3*a+b)/4
          test = if bflag then abs(s - b) < (1/2)*abs(b - c) else c == d || abs(s - b) < (1/2)*abs(c - d)

--Printing wrapper for Results
printClean :: String -> Result -> IO()
printClean method res = putStrLn $ method ++ " took " ++ show(steps res) ++ " steps to find "  ++ show(realToFrac(root res))

main = do
    let a = 1
    let b = 2
    printClean "Bisection Method" (bisect a b)
    printClean "Secant Method" (secant b a)
    printClean "Dekker's Method" (dekker a b)
    printClean "IQI" (invquadinterp (realToFrac b) (realToFrac ((a+b)/2)) (realToFrac a))
    printClean "Brent's Method" (brentdekkerill a b)
