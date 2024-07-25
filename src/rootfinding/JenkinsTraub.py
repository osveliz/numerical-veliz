from sympy import Poly, symbols, div, I, nan, zoo
import numpy as np
import matplotlib.pyplot as plt

x = symbols('x')

"""
Implementation of Traub's and simplified Jenkins-Traub algorithms

Author
------
Oscar Veliz
Date: 2024-07-25
Version: 1.0
"""


def compute_G(P: Poly, n=10):
    """
    Computes Traub's largest root finder using G polynomial.

    Parameters
    ----------
    P : Poly
        The input polynomial.
    n : int, optional
        The number of iterations to perform (default is 10).

    Returns
    -------
    Poly
        The polynomial G after n iterations
    """
    G = [P.diff(x)]
    for _ in range(n):
        lc = G[-1].LC()
        G_next = x * G[-1] - lc * P
        G.append(G_next)
    # print(G) # uncomment to see all G
    return G[-1]


def shift(P: Poly, H: Poly, s=0):
    """
    Shifts a polynomial H by evaluating it at a shift point s using polynomial division.

    Parameters
    ----------
    P : Poly
        The input polynomial.
    H : Poly
        The polynomial to be shifted.
    s : float, optional
        The shift point (default is 0).

    Returns
    -------
    Poly
        The shifted H.
    """
    Ps = P.eval(s).evalf()
    Hs = H.eval(s).evalf()
    result, _ = div(H-P*Hs/Ps, x-s)
    return result


def compute_H(P: Poly, n=10, s=0):
    """
    Computes Traub's H polynomial.

    Parameters
    ----------
    P : Poly
        The input polynomial.
    n : int, optional
        The number of iterations to perform (default is 10).
    s : float, optional
        The shift point for evaluating H (default is 0).

    Returns
    -------
    Poly
        The polynomial H after n iterations.
    """
    H = [P.diff(x)]
    for _ in range(n):
        H.append(shift(P, H[-1], s))
    # print(H) # uncomment to see all H
    return H[-1]


def runner(P: Poly, GH: Poly, s=1.1, eps=10**-7):
    """
    Runs Traub's method for given P and either a G or H poly at a given point s.

    This function replaces normalized G or H instead of derivative in Newton's Method.

    Parameters
    ----------
    P : Poly
        The input polynomial.
    GH : Poly
        G or H polynomial.
    s : float, optional
        The starting point for root finding (default is 1.1).
    eps : float, optional
        The convergence tolerance (default is 10^-7).

    Returns
    -------
    float
        The approximated root of the polynomial.
    """
    Px = P.eval(s).evalf()
    GHx = GH.eval(s).evalf()
    lc = GH.LC()
    i = 0
    while abs(Px) > eps and i < 50 and GHx != nan and GHx != zoo and Px != nan and Px != zoo:
        s = s - Px / (GHx / lc)
        Px = P.eval(s).evalf()
        GHx = GH.eval(s).evalf()
        i += 1
        # print("s",s) # uncomment to see iterations
    if i == 50:
        return nan
    return s.evalf()


def traub1(P: Poly, x=1, L=10, eps=10**-7):
    """
    Finds the largest modular root of a polynomial using Traub's G polynomial.

    Parameters
    ----------
    P : Poly
        The input polynomial.
    x : float, optional
        The starting point for root finding (default is 1).
    L : int, optional
        The number of iterations for computing G (default is 10).
    eps : float, optional
        The convergence tolerance (default is 10^-7).

    Returns
    -------
    float
        The approximated largest root of the polynomial.
    """
    G = compute_G(P, L)
    return runner(P, G, x, eps)


def traub2(P: Poly, x=0, L=10, eps=10**-7):
    """
    Finds the smallest modular root of a polynomial using Traub's H polynomial.

    Parameters
    ----------
    P : Poly
        The input polynomial.
    x : float, optional
        The starting point for root finding (default is 0). To see what JenkinsTraub
        would do with a fixed shift stage 2, replace the default x value.
    L : int, optional
        The number of iterations for computing H (default is 10).
    eps : float, optional
        The convergence tolerance (default is 10^-7).

    Returns
    -------
    float
        The approximated smallest root of the polynomial.
    """
    H = compute_H(P, n=L)
    return runner(P, H, x, eps)


def stage2(P: Poly, H: Poly, s: float, M=5):
    """
    Performs Stage 2 of the Jenkins-Traub method.

    Parameters
    ----------
    P : Poly
        The input polynomial.
    H : Poly
        The polynomial from Stage 1.
    s : float
        The shift value s.
    M : int, optional
        The number of iterations (default is 5).

    Returns
    -------
    Poly
        The updated H polynomial.
    """
    H = [H]
    for _ in range(M):
        Hx = H[-1]
        Hs = Hx.eval(s).evalf()
        if Hs == nan or Hs == zoo:
            H.pop()
            break
        H_next = shift(P, H[-1], s)
        H.append(H_next)
    # print("stage 2 H", H)  # uncomment to see all H
    return H[-1]


def stage3(P: Poly, H: Poly, s: float, withCount=False, eps=10**-7):
    """
    Performs Stage 3 of the Jenkins-Traub method.

    Parameters
    ----------
    P : Poly
        The input polynomial.
    H : Poly
        The polynomial from Stage 2.
    s : float
        The initial guess.
    withCount : bool, optional
        If True, returns the number of iterations instead of the root (default is False).
    eps : float, optional
        The convergence tolerance (default is 10^-7).

    Returns
    -------
    float
        The final approximated root of the polynomial, or the number of iterations if withCount is True.
    """
    H = [H]
    Ps = P.eval(s).evalf()
    HBar = normH(H[-1])
    HBs = HBar.eval(s).evalf()
    s = s - Ps/HBs
    i = 0
    while abs(Ps) > eps and i < 25:
        H_next = shift(P, H[-1], s)
        HBar = normH(H[-1])
        HBs = HBar.eval(s).evalf()
        s = s - Ps/HBs
        Ps = P.eval(s).evalf()
        H.append(H_next)
        i += 1
    # print("stage 3 H", H) # uncomment to see all H
    return i if withCount else s.evalf()


def normH(p: Poly):
    """
    Normalizes the polynomial p by dividing by its leading coefficient.

    Parameters
    ----------
    p : Poly
        The polynomial to be normalized.

    Returns
    -------
    Poly
        The normalized polynomial.
    """
    return 1/p.LC() * p  # if used p/p.LC() it would need to be cast as a Poly


def simpleJK(P: Poly, s=1.1, withCount=False):
    """
    Performs a simplified Jenkins-Traub method to find a root of a polynomial.

    This method includes a basic Stage 1 but not including the convoluted way
    of approximating a good s value via Cauchy's beta method.

    Parameters
    ----------
    P : Poly
        The input polynomial.
    s : float, optional
        The starting point for root finding (default is 1.1).
    withCount : bool, optional
        If True, returns the number of iterations instead of the root (default is False).

    Returns
    -------
    float
        The approximated root of the polynomial, or the number of iterations if withCount is True.
    """
    H = [P.diff(x)]  # basic stage 1
    # H = [compute_H(P,3)] # towards smallest stage 1
    # print("stage 1 H", H[-1])  # uncomment to see stage 1
    Hx = stage2(P, H[-1], s, 1)
    # print("stage 2 H", Hx)  # uncomment to see stage 2
    return stage3(P, Hx, s, withCount)


def fractal(P: Poly, xmin=-3, xmax=3, ymin=-1.69, ymax=1.69):
    """
    Create a fractal based on Jenkins-Traub.

    The horizontal represtents real numbers while vertical are imaginary.
    P is a polynomial. Saved to "JenkinsTraub.png".

    Parameters
    ----------
    P : Poly
        The polynomial coefficients p(x) = ax² + bx + c --> [a,b,c].
    xmin : float, optional
        Leftmost number. Default -3.
    xmax : float, optional
        Rightmost number. Default 3.
    ymin : float, optional
        Bottommost number. Default -1.69.
    ymax : float, optional
        Topmost number. Default 1.69.
    """
    plt.figure(frameon=False)
    fig, ax = plt.subplots()
    fig.set_size_inches(19.1, 10.75)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    real, img = np.meshgrid(np.linspace(ymin, ymax, 2160),  # 3840 x 2160 is very hi-res takes a while
                            np.linspace(xmin, xmax, 3840))
    vjk = np.vectorize(simpleJK, excluded=['a'])
    z = vjk(P, real + img * I, withCount=True)
    ax.pcolormesh(real, img, z, vmin=z.min(), vmax=z.max(),
                  cmap=plt.get_cmap('ocean'), shading='nearest')
    plt.savefig("JenkinsTraub.png", transparent=True)


def main():
    """
    Execute the Traub and simplified Jenkins-Traub methods.
    """
    # P = Poly(x**2 - x - 1)
    P = Poly(x**5 - 8*x**4 - 72*x**3 + 382*x**2 + 727*x - 2310)
    # P = Poly(x**3 - 1.0)
    s = 4.1
    print(P)
    print("Traub G", traub1(P, s))
    print("Traub H", traub2(P, s))
    print("Jenkins", simpleJK(P, s))
    # fractal(P)  # uncomment to see fractal
    print("done")


if __name__ == "__main__":
    main()
