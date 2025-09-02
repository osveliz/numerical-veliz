"""
Lehmer-Schur Implementation for the Oscar Veliz YouTube channel

Uses symbolic programming to make understanding the functions easier. This program
also draws circles in the complex plane visualizing the search. Simply change the
function string in main to be whatever polynomial you wish.

Author: Oscar Veliz

References:
    [1] Lehmer, D.H. "A machine method for solving polynomial equations." Journal of the ACM 8.2 (1961).
    [2] Lehmer, D.H. "Search procedures for polynomial equation solving." Constructive Aspects of the Fundamental Theorem of Algebra (1969). 
"""

from sympy import symbols, Poly, conjugate, exp, pi, sympify, latex
import matplotlib.pyplot as plt

x = symbols('x')
fig, ax = plt.subplots(figsize=(8, 8))


def draw_circle(c: complex, r: float, found: bool):
    circle = plt.Circle((c.real, c.imag), r, linewidth=0.5,
                        facecolor='none', edgecolor='black',
                        linestyle='-' if found else '--')
    ax.add_patch(circle)


def calculate_g(f: Poly, r: float, c: complex):
    g = Poly(f.subs(x, r * x + c), x)
    return g.simplify()


def calculate_gstar(g: Poly):
    a = g.all_coeffs()
    new_a = [conjugate(coeff) for coeff in a]
    g_star = sum(coef * x**i for i, coef in enumerate(new_a))
    return g_star.simplify()


def calculate_T(g: Poly):
    g_star = calculate_gstar(g)
    an = g.LC()
    a0 = g.coeff_monomial(x**0)
    T = Poly(conjugate(a0)*g - an*g_star)
    return T.simplify()


def contains_root(f: Poly, r: float, c: complex):
    g = calculate_g(f, r, c)
    T = calculate_T(g)
    g0 = T.eval(0).evalf()
    while g0 > 0:
        T = calculate_T(T)
        g0 = T.eval(0).evalf()
    return g0 < 0


def search(f: Poly, r=1.0, c=0+0j, eps=10**-4):
    c = complex(c)  # just making sure c is complex
    if abs(f.eval(c).evalf()) < eps:
        plt.plot(c.real, c.imag, 'k.')
        return c
    if contains_root(f, r, c):  # start has root - shrink
        upper_r = r
        draw_circle(c, upper_r, True)
        lower_r = upper_r/2
        while contains_root(f, lower_r, c):
            upper_r = lower_r
            draw_circle(c, upper_r, True)
            lower_r = lower_r/2
        draw_circle(c, lower_r, False)
    else:  # start does not have root - expand
        lower_r = r
        draw_circle(c, lower_r, False)
        upper_r = lower_r*2
        while not contains_root(f, upper_r, c):
            lower_r = upper_r
            draw_circle(c, lower_r, False)
            upper_r = upper_r * 2
        draw_circle(c, upper_r, True)
    cs = [5 * lower_r / 3 * exp(2j * pi * (i * 3 %
                                8) / 8).evalf() + c for i in range(8)]
    found = False
    i = 0
    shrunk_r = 5*lower_r/6
    while not found:
        currentc = cs[i]
        found = contains_root(f, shrunk_r, currentc)
        draw_circle(complex(currentc), shrunk_r, found)
        i = i + 1
    return search(f, shrunk_r/2, currentc)


if __name__ == '__main__':
    p = '16*x**4-8*x**2+9'
    f = sympify(p).as_poly()
    print('f(x): ', f)
    root = search(f)
    print('r = ', root, 'f(r) = ', f.eval(root).evalf())

    ax.set_aspect('equal')
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3, 3)
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.title('Lehmer-Schur $'+latex(f.as_expr())+'$')
    plt.xlabel('Real')
    plt.ylabel('Imaginary')
    plt.show()
