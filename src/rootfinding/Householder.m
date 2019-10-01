#################################
# GNU Octave program comparing 
# Newton, Halley, and Householder
# created for the YouTube channel
# https://youtube.com/OscarVeliz
# @author Oscar Veliz
#################################
pkg load symbolic

################################
# Householder's Method
# @param f function string
# @param xn starting value
# @param d derivative degree
# return root of f
#################################
function root = householder(f, xn, d)
    syms x;
    maximum = 10;
    epsilon = 10 ^ -12;
    fx = function_handle(sym(f));
    # check if continuously differentiable
    continuous = char(simplify(sym(f)));
    for i = 1:d
        continuous = char(simplify(diff(sym(continuous))));
        if strcmp(continuous,"0")
            printf("f only %d times continuously differentiable \n",i-1)
            d = i-1;
            break;
        endif
    endfor
    # find (1/f)^(d) derivatives
    derivs = cell(1,d);
    oneoverf = strcat("1.0 / (" , f , ")");
    derivs(1) = char(simplify(diff(sym(oneoverf))));
    for i = 2:d
        derivs(i) = char(simplify(diff(sym(derivs(i-1)))));
    endfor
    derivs; # remove ; to see derivatives
    # create functions from derivatives
    if (d == 1)
        numer = sym(oneoverf);
    else
        numer = sym(derivs(d-1));
    endif
    denom = sym(derivs(d));
    eq = simplify(sym(strcat("(",char(numer),")/(",char(denom),")")))
    house = function_handle(eq);
    # perform iteration
    printf("x0 = %.20g\n",xn*1.0)
    c = 0;
    while and(abs(fx(xn)) > epsilon , (c < maximum))
        xn = xn + d * house(xn);
        c++;
        printf("x%d = %.20g\n",c,xn)
    endwhile
    root = xn;
endfunction

################################
# Newton's Method
# @param f function string
# @param xn starting value
# return root of f
#################################
function root = newton(f, xn)
    root = householder(f, xn, 1);
endfunction

################################
# Halley's Method
# @param f function string
# @param xn starting value
# return root of f
#################################
function root = halley(f, xn)
    root = householder(f, xn, 2);
endfunction

# main
func = "x^8 + 15*x^4 - 16"
begin = 2
printf("\nNewton's Method\n")
newton(func, begin)
printf("\nHalley's Method\n")
halley(func, begin)
printf("\nHouseholder's Method - 3\n")
householder(func, begin, 3)
printf("\nHouseholders's Method\n")
householder(func, begin, 10)
