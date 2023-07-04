import math  # to compare with built in e

"""
An implementation of e-spigot algorithm
by Rabinowitz & Wagon 1995 for computing
Euler's Number to any degree of accuracy

@author Oscar Veliz
"""

n = 1006  # number of desired digits (add 6 for safety)
A = n*[1]
e = "2."  # start with leading 2
for s in range(n):
    for i in range(n):
        A[i] = A[i] * 10
    for i in range(n):
        j = n-i-1  # start from right side
        k = j + 2  # leftmost division is by 2 not 0
        r = A[j] % k  # remainder
        q = A[j] // k  # quotient
        if j != 0:
            A[j-1] = A[j-1] + q  # carry the quotient
        else:
            e = e + str(q)  # last quotient is digit of e
        A[j] = r
print(math.e) # built-in to compare
print(e[:-6]) # remove extra 6 spots
