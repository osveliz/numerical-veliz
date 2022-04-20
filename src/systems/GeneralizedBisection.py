"""Generalized Bisection Visualization

Method of Bisecting Triangles to solve a System
of Nonlinear Equations using simplified version
of Harvey-Stenger 2D Analogue Method.

Requires numpy and matplotlib.

:author: Oscar Veliz
"""
import numpy as np
import matplotlib.pyplot as plt


def F(x, y):
    """Function of 2 variables
    :param x: the x value
    :param y: the y value
    :return: tuple. (x²-y-1,x-y²+1)
    """
    return (x**2.0-y-1.0, x-y**2.0+1.0)


def npF(X):
    """numpy version of F"""
    return np.array(F(X[0], X[1]))


def triToPoints(T):
    """Convert Matrix of 3 points to tuple
    :param T: Triangle with 3 points
    :return: (A,B,C)
    """
    return (T[0, :], T[1, :], T[2, :])


def L(A, B, X):
    """Harvey-Stenger linear form L(A,B,X)
    :param A: [x,y] point as ndarray
    :param B: [x,y] point as ndarray
    :param X: [x,y] point as ndarray
    :return: result of linear form equation
    """
    return (B[1]-A[1])*(X[0]-A[0]) - (B[0]-A[0])*(X[1]-A[1])


def check(T, V):
    """Harvey-Stenger test if a point V is inside a triangle with points in the matrix T
    :param T: Three points A, B, and C in a matrix
    :param V: Point to test
    :return: True when V is in T, False otherwise
    """
    A, B, C = triToPoints(T)
    return L(A, B, V)*L(A, B, C) >= 0 and L(B, C, V)*L(B, C, A) >= 0 and L(C, A, V)*L(C, A, B) >= 0


def center(T):
    """Compute the center of a Triangle
    :param T: A triangle
    :return: the average (center) of the three points
    """
    A, B, C = triToPoints(T)
    return (A+B+C)/3.0


def eval(T):
    """Evaluate a matrix of 3 points
    :param T: Matrix of 3 points A, B, and C
    :return: [F(A),F(B),F(C)] as ndarray
    """
    A, B, C = triToPoints(T)
    return np.array([npF(A), npF(B), npF(C)])


def rotate(T):
    """Given matrix of three points, rotate the points so AB is longest
    :param T: Matrix representing a triangle with 3 points
    :return: Rotated version of same Triangle
    """
    A, B, C = triToPoints(T)
    ab = np.linalg.norm(A-B)
    bc = np.linalg.norm(B-C)
    ca = np.linalg.norm(C-A)
    if ab >= bc and bc >= ca:
        return np.array([A, B, C])
    elif ab >= ca and ca >= bc:
        return np.array([B, A, C])
    elif bc >= ca and ca >= ab:
        return np.array([B, C, A])
    elif bc >= ab and ab >= ca:
        return np.array([C, B, A])
    elif ca >= ab and ab >= bc:
        return np.array([C, A, B])
    else:
        return np.array([A, C, B])


def drawTri(T):
    """Draws a triangle given matrix with three points"""
    s = np.vstack([T, T[0, :]])
    x = s[:, 0]
    y = s[:, 1]
    plt.plot(x, y, '-')
    plt.fill(x, y, alpha=0.1)


def setup(d=(-3.5, 3.5), r=(-2, 2), size=(16, 9), res=100000):
    """Sets up the plotting space
    :param d: domain tuple default to (-3.5,3.5)
    :param r: range tuple default to (-2,2)
    :param size: figure size default to (16,9)
    :param res: number of points in plot, default 1000000
                Lower the res if it is taking too long to plot
    """
    x = np.linspace(d[0], d[1], res)
    y = x**2 - 1  # x² - y - 1
    plt.figure(figsize=size)
    plt.plot(x, y, 'b', label='x²-y-1')
    plt.plot(y, x, 'g', label='x-y²+1')
    plt.grid(True, linestyle=':')
    plt.xlim([d[0], d[1]])
    plt.ylim([r[0], r[1]])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.title('Generalized Bisection')


def main():
    """Implementation of 2D Bisection based on simplified
    version of Harvey-Stenger 2D Analogue
    """
    setup()
    A = np.array([1.0, 0.5])
    B = np.array([1.5, 2.0])
    C = np.array([2.0, 2.0])
    D = np.array([2.0, 1.5])
    R = np.array([A, B, C])
    S = np.array([A, D, C])
    drawTri(R)
    drawTri(S)
    i = 0
    eps = 10**-6
    zero = np.zeros(2)
    E = (A+B+C+D)/4
    longest = max(np.linalg.norm(C-A), np.linalg.norm(D-B))
    while np.linalg.norm(npF(E)) >= eps and longest >= eps and i < 100:
        evalR = eval(R)
        evalS = eval(S)
        # drawTri(evalR) # uncomment to see searched area
        # drawTri(evalS) # uncomment to see searched area
        if check(evalR, zero):
            T = R
            # print("R") # uncomment to see choice
        elif check(evalS, zero):  # change to else for efficiency
            T = S
            # print("S") # uncomment to see choice
        else:  # for safety, when neither contains zero
            print("problem")
            break
        T = rotate(T)
        drawTri(T)
        A, B, C = triToPoints(T)
        E = center(T)
        D = (A+B)/2.0
        longest = np.linalg.norm(B-A)
        R = np.array([A, D, C])
        S = np.array([D, B, C])
        i += 1
    print(E, np.linalg.norm(npF(E)), i)
    print(T)
    plt.plot(0, 0, marker="+", markersize=7, markeredgecolor="black")
    #plt.plot(E[0], E[1], marker=".",markersize=7,markerfacecolor="black", markeredgecolor="black")
    plt.show()


if __name__ == "__main__":
    main()
