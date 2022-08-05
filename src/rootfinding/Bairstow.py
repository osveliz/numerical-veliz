import numpy as np
import matplotlib.pyplot as plt
"""Bairstow's Method Implementation.

Polynomial Solver, including complex roots, by extracting a quadratic
factor (x² - ux - v) from a quotient until the remainder is zero. Then
repeat with the quotient as the new polymonial until left with degree
2 or less then solve directly. See lesson by Oscar Veliz on YouTube.

@author: Oscar Veliz

Notes
-----
.. [1] Bairstow, Leonard. Applied aerodynamics. Longmans, Green and Company, 1920.
   https://books.google.com/books?id=GIUQAQAAMAAJ
.. [2] Henrici, Peter. "Elements of numerical analysis." (1964).
   https://archive.org/details/elements-of-numerical-analysis-by-peter-henrici
"""

def horner(p, x):
    """Horner-Ruffini Method for evaluating polynomials.
    
    The list p contains the coefficients [a, b, c, ...]:

    p(x) = ax² + bx + c --> ((a)x + b)x + c

    This function is only used to check that the found roots result
    in near zero values. It isn't used as part of Bairstow's Method,
    rather it is a fast way to evaluate p(x).

    Parameters
    ----------
    p : list of float
        The polynomial coefficients p(x) = ax² + bx + c --> [a,b,c].
    x : float or complex
        The number to evaluate i.e. p(3) = ?.

    Returns
    -------
    float
        The result of p(x).
    """
    result = p[0]
    for i in range(1, len(p)):
        result = p[i] + x*result
    return result


def quo(a, u, v):
    """Determine the quotient of a polynomial.
    
    Coefficients given in list a and terms u and v of dividing quadtratic:

    p(x) = a_0*x^n + a_1*x^{n-1} + ... + a_{n-1}x^1 + a_n

    q(x) = (b_0*x^{n-2} + b_1*x^{n-3} + ... + b_{n-2})(x²-ux-v) + b_{n-1}x + b_n

    When b_{n-1} and b_n are both zero then u and v make a quadratic
    with roots of p(x) since there is no remainder.

    Parameters
    ----------
    a : list of float
        The polynomial coefficients p(x) = ax² + bx + c --> [a,b,c].
    u : float
        The coefficient in (x²-ux-v).
    v : float
        The coefficient in (x²-ux-v).

    Returns
    -------
    list of float
        Quotient q(x) as list of coefficients including remainder.
    """
    b = [0, 0]
    for i in range(len(a)):
        b.append(a[i] + u*b[i+1] + v*b[i])
    return b[2:]


def quadratic(u, v):
    """Quadratic equation for x² - ux - v for finding roots.

    Note that the usual quadratic uses positive signs ax² + bx + c where:
    
    roots = (-b ± sqrt(b² - 4ac))/(2a)

    This differs from this function since a is 1 and both u and v have
    opposite signs. This makes our form:

    roots = (u ± sqrt(u² + 4v))/2

    Parameters
    ----------
    u : float
        The coefficient in x² - ux - v.
    v : float
        The coefficient in x² - ux - v.

    Returns
    -------
    tuple of two floats.
        Roots of given polynomial u=0 v=1 --> roots = (1,-1).
    """
    disc = (u**2+4*v)**(0.5)
    return ((u+disc)/2, (u-disc)/2)


def bairstow(a, u, v, eps=10**(-12), max=50, count=False):
    """Find quotient of polynomial using Bairstow's method.
    
    Coefficients given in list a and terms u and v of quadtratic.
    It iterates on u and v until it finds a pair that cause the
    remainder to be zero.

    p(x) = a_0*x^n + a_1*x^{n-1} + ... + a_{n-1}x^1 + a_n

    p(x) = q(x)(x²-ux-v) + remainder

    When remainder is zero (last two terms in result from quo) are both zero,
    stop iterating. Final results are the coefficients from quo and final u and v.

    Parameters
    ----------
    a : list of float
        The polynomial coefficients p(x) = ax² + bx + c --> [a,b,c].
    u : float
        Initial value for the coefficient in (x²-ux-v).
    v : float
        Initial value for the coefficient in (x²-ux-v).
    eps : float, optional
        Default epsilon value is 10^{-12}.
    max : int, optional
        Default maximum iteration count is 50.
    count : bool, optional
        Default to false. When True, return the number of iterations
        it took to find final values for u and v.

    Returns
    -------
    tuple or integer
        When count is True, return the number of iterations to find solution.
        Otherwise return a tuple containing new coefficients of quotient and
        values for u and v.
    
    See Also
    --------
    quo : Quotient.
    """
    i = 0
    n = len(a) - 1  # length of c
    b = quo(a, u, v)
    while (abs(b[-1]) > eps or abs(b[-2]) > eps) and i < max:
        c = quo(b, u, v)[:-1]
        denom = c[n-2]**2 - c[n-1]*c[n-3]
        du = (b[n]*c[n-3] - b[n-1]*c[n-2]) / denom  # b_n-1 typo in henrici
        dv = (b[n-1]*c[n-1] - b[n]*c[n-2]) / denom
        u += du
        v += dv
        i += 1
        b = quo(a, u, v)
    return i if count else  (b[:-2], u, v)

def allRoots(p):
    """Find all roots of a polynomial using Bairstow's Method.

    Once a quotient is found, solve the quadratic and add to list of roots.
    When degree of polynomial is 2 or less, solve directly and add to list.

    Parameters
    ----------
    p : list of float
        The polynomial coefficients p(x) = ax² + bx + c --> [a,b,c].

    Returns
    -------
    list of float
        All roots of p(x).
    """
    a = [p[i]/p[0] for i in range(len(p))]  # normalize
    print("a =", a)
    roots = []
    while(len(a) > 3):
        u = 1
        v = 1
        a, u, v = bairstow(a, u, v)
        print(a)
        print("u =",u,"v =",v)
        r1, r2 = quadratic(u, v)
        print("r =",r1,r2)
        roots.append(r1)
        roots.append(r2)
    if(len(a) == 3):  # degree 2
        print(a)
        print("u =",-a[1],"v =",-a[2])
        r1, r2 = quadratic(-a[1], -a[2])
        print("r =",r1,r2)
        roots.append(r1)
        roots.append(r2)
    else: # degree 1
        print(a)
        print("r =",-a[1])
        roots.append(-a[1])
    print("r =", roots)
    print("p(r) =", [horner(p, r) for r in roots])
    return roots


def fractal(p, xmin=-3, xmax=3, ymin=-2, ymax=2,):
    """Create a fractal based on Bairstow's Method.
    
    The horizontalreprestents values for u, while vertical are v values.
    The list p contains the coefficients. Saved to "BairstowFractal.png".

    Parameters
    ----------
    p : list of float
        The polynomial coefficients p(x) = ax² + bx + c --> [a,b,c].
    umin : float, optional
        Leftmost number. Default -3.
    umax : float, optional
        Rightmost number. Default 3.
    vmin : float, optional
        Bottommost number. Default -2.
    vmax : float, optional
        Topmost number. Default 2.
    """
    plt.figure(frameon=False)
    fig, ax = plt.subplots()
    fig.set_size_inches(16, 10.75)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    v, u = np.meshgrid(np.linspace(ymin, ymax, 1000), np.linspace(xmin, xmax, 2000))
    vbairstow = np.vectorize(bairstow, excluded=['a'])
    z = vbairstow(a=p, u=u, v=v, count=True)
    ax.pcolormesh(u, v, z, vmin=z.min(), vmax=z.max(), cmap=plt.get_cmap('ocean'), shading='nearest')
    plt.savefig("BairstowFractal.png", transparent=True)


def main():
    """Main Function."""
    allRoots([1, 0, 0, 0, 15, 0, 0, 0, -16])
    allRoots([1, 0, 0, -1])
    allRoots([1, -8, -72, 382, 727, -2310])
    allRoots([1, -1, 2, 5])
    allRoots([1,20.4,151.3,490,687,719,150,109,6.87])
    fractal([1, 0, 0, 0, 1])


if __name__ == "__main__":
    main()
