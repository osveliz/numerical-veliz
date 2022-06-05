#################################
# GNU Octave Symbolic Program
# Graeffe's Root-Squaring Method
# Created for the YouTube channel
# https://youtube.com/OscarVeliz
# @author Oscar Veliz
#################################
pkg load symbolic

################################
# Create a new polynomial whose
# roots are squared from original
# @param f function string
# @return a new polynomial with squared roots
#################################
function newpoly = enke(f)
    syms x y;
    [c,t] = coeffs(f);
    right = sym();
    left = sym();
    for i = 1:length(c)
        deg = degree(t(i));
        if mod(deg, 2) == 0
            right = right - c(1, i)*x^deg;
        else
            left = left + c(1, i)*x^deg;
        endif
    endfor
    left = left^2;
    right = right^2;
    newpoly = expand(left - right);
    newpoly = subs(newpoly, x^2, y);
    newpoly = subs(newpoly, y, x);
endfunction

################################
# Approximate the roots of a function given m
# @param f function string
# @param m exponent of roots
# @return a new polynomial with squared roots
#################################
function guess = approx(f, m)
    [c,t] = coeffs(f);
    n = length(c);
    deg = degree(f);
    if n != degree(f) + 1 # missing term(s)
        guess = ones(1, deg);
        return;
    endif
    logs = sym([]);
    guess = sym([]);
    logs(1) = log10(abs(c(2)))/m;
    guess(1) = vpa(10^logs(1));
    for i = 2:n-1;
        logs(i) = log10(abs(c(i + 1)))/m;
        guess(i) = vpa(10^(logs(i) - logs(i - 1)));
    endfor;
endfunction

################################
# Creates Wilkinson Polynomial
# @return Wilkinson Polynomial
#################################
function w = wilkinson()
    syms x;
    w = (x - 1);
    for i = 2:20
        w = w*(x - i);
    endfor
endfunction

################################
# Guesses are in absolute values,
# find if the guess should be +/-
# @param guess - array of guesses
# @return array containing roots
#################################
function solutions = testGuesses(f, guess)
    solutions = sym([]);
    fx = function_handle(sym(f));
    for i = 1:length(guess)
        g = guess(i);
        if abs(fx(g)) < abs(fx(-g))
            solutions(i) = g;
        else
            solutions(i) = -g;
        endif
    endfor
endfunction

################################
# Graeffe - Dandelin - Lobachevskiĭ
# @return array of absolute valued roots
#################################
function guess = gdl(f)
    lastguess = sym(zeros(1, degree(f)));
    guess = approx(f, 1)
    i = 0;
    k = 10;
    eps = sym(1)/10^6;
    while (i < k && norm(guess - lastguess) > eps)
        i = i + 1
        f = enke(f)
        lastguess = guess;
        guess = approx(f, 2^i)
    endwhile
endfunction

### Main
syms x;
f = x^3 - 13*x^2 + 32*x - 20
[c,t] = coeffs(f)
absroots = gdl(f)
solutions = testGuesses(f, absroots)
#answers = roots(sym2poly(f)) # uncomment to see answers
