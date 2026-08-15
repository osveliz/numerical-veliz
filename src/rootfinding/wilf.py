"""
Wilf's Method aka Global Bisection uses
Sturm-sequences for finding all roots of a
polynomial including complex roots.

Created for the Oscar Veliz YouTube channel
"""

import sympy as sp
import matplotlib.pyplot as plt

x = sp.symbols('x')
t = sp.symbols('t', real=True)
z = sp.symbols('z')

EPS = 1e-9
CMAP = plt.get_cmap('Paired')


def starting_square(coeffs):
    """Initial square 'radius' from Cauchy.
    Center offset by small amount so bisections never land exactly on 
    the real or imaginary axis."""
    R = 1 + max(abs(c / coeffs[0]) for c in coeffs[1:])  # Cauchy Bound
    x0, x1 = -R + 1e-5, R + 1.3e-5
    y0, y1 = x0, x1
    return [complex(x1, y1), complex(x0, y1),
            complex(x0, y0), complex(x1, y0)]


def sign(x):
    if x > 1e-11:
        return 1
    if x < -1e-11:
        return -1
    return 0


def V(sequence, point):
    """Sign-variation count of the sequence."""
    vals = []
    for poly in sequence:
        s = sign(float(sp.N(poly.as_expr().subs(t, point))))
        if s != 0:
            vals.append(s)
    return sum(1 for i in range(len(vals) - 1) if vals[i] != vals[i + 1])


def edge_interval_calculation(P_expr, Q_start, Q_end):
    """V(L) - V(0) for one side of the square."""
    Q_start, Q_end = sp.sympify(Q_start), sp.sympify(Q_end)
    L = sp.sqrt((sp.re(Q_end) - sp.re(Q_start)) ** 2 +
                (sp.im(Q_end) - sp.im(Q_start)) ** 2)
    if L == 0:
        return 0
    direction = (Q_end - Q_start) / L
    z_t = Q_start + t * direction
    poly_t = sp.Poly(sp.expand(P_expr.subs(z, z_t)), t)
    coeffs_complex = poly_t.all_coeffs()
    f1 = sp.expand(sum(sp.re(c) * t ** i for i,
                   c in enumerate(reversed(coeffs_complex))))
    f2 = sp.expand(sum(sp.im(c) * t ** i for i,
                   c in enumerate(reversed(coeffs_complex))))
    seq = [sp.Poly(f1, t), sp.Poly(f2, t)]
    while seq[-1].degree() > 0:
        prev2, prev1 = seq[-2], seq[-1]
        remainder = sp.rem(prev2.as_expr(), prev1.as_expr(), t)
        nxt = sp.Poly(sp.expand(-remainder), t)
        seq.append(nxt)
        if nxt.is_zero:
            break
    return V(seq, L) - V(seq, 0)


def wilfs_formula(P_expr, corners):
    """N(P; R)."""
    n = len(corners)
    return sum(edge_interval_calculation(P_expr, corners[k], corners[(k + 1) % n])
               for k in range(n)) // 2


def color_for_n(n, max_n):
    if n == 0:
        return (0.6, 0.6, 0.6, 0.5)  # light gray, empty box
    return CMAP(1.0 if max_n <= 1 else (n - 1) / (max_n - 1))


def draw_rectangle(ax, x0, x1, y0, y1, n, max_n):
    lw = 0.5 if n == 0 else 1.5 + 0.5 * n
    rect = plt.Rectangle((x0, y0), x1 - x0, y1 - y0, linewidth=lw,
                         facecolor='none', edgecolor=color_for_n(n, max_n),
                         linestyle='--' if n == 0 else '-')
    ax.add_patch(rect)


def bisect(ax, P_expr, x0, x1, y0, y1, max_n, max_depth, depth=0):
    """Recursively quarters the search space; stop when squares become small enough."""
    n = wilfs_formula(P_expr, [complex(x1, y1), complex(x0, y1),
                               complex(x0, y0), complex(x1, y0)])
    draw_rectangle(ax, x0, x1, y0, y1, n, max_n)
    if n == 0:
        return []
    if max(x1 - x0, y1 - y0) < EPS or depth > max_depth:
        c = complex((x0 + x1) / 2, (y0 + y1) / 2)
        ax.plot(c.real, c.imag, 'k.')
        return [c]
    xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
    roots = []
    for qx0, qx1, qy0, qy1 in [(x0, xm, y0, ym), (xm, x1, y0, ym),
                               (x0, xm, ym, y1), (xm, x1, ym, y1)]:
        roots += bisect(ax, P_expr, qx0, qx1, qy0, qy1, max_n,
                        max_depth, depth + 1)
    return roots


def search(ax, coeffs, max_depth):
    """Set up for Wilf's Method and plotting"""
    P_expr = sum(c * z ** i for i, c in enumerate(reversed(coeffs)))
    max_n = len(coeffs) - 1
    corners = starting_square(coeffs)
    xs = [c.real for c in corners]
    ys = [c.imag for c in corners]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    roots = bisect(ax, P_expr, x0, x1, y0, y1, max_n, max_depth)
    return roots, (x0, x1, y0, y1), max_n


def main():
    coeffs = [1, 0, 0, -1]  # x^3 - 1
    # coeffs = [16, 0, -8, 0, 9]  # 16x^4 - 8x^2 + 9
    print('coeffs:', coeffs)
    fig, ax = plt.subplots(figsize=(8, 8))
    roots, (x0, x1, y0, y1), max_n = search(ax, coeffs, 50)
    print(f'Found {len(roots)} approximate root locations:')
    for r in roots:
        print(f'  {r}')

    ax.set_aspect('equal')
    pad = 0.05 * (x1 - x0)
    ax.set_xlim(x0 - pad, x1 + pad)
    ax.set_ylim(y0 - pad, y1 + pad)
    ax.set_title(r'Wilf $' + sp.latex(sum(c * x ** i for i,
                 c in enumerate(coeffs[::-1]))) + '$')
    ax.set_xlabel('Real')
    ax.set_ylabel('Imaginary')
    sm = plt.cm.ScalarMappable(
        cmap=CMAP, norm=plt.Normalize(vmin=1, vmax=max(max_n, 1)))
    sm.set_array([])
    fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.04, ticks=range(1, max_n + 1))
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
