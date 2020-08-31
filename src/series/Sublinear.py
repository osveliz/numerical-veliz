###########################
# Sublinear Convergence of
# Gregory-Leibniz Series
# See Machin.lsp
# @author Oscar Veliz
###########################
import math
itr = 30000
series = 0.0
terms = []
e = []
for i in range(itr):
    series += pow(-1,i) / (2*i+1)
    terms.append(series)
    e.append(abs(4*series - math.pi))
a = [0]
m = [0]
print("n ~π e a m")
print(0,4*terms[0],e[0])
for i in range(1,itr-1):
    a.append(math.log(e[i+1]/e[i]) / math.log(e[i]/e[i-1]))
    m.append(e[i+1]/(math.pow(e[i],a[i])))
    print(i,4*terms[i],e[i],a[i],m[i])
print(itr-1,4*terms[itr-1],e[itr-1],'-','-')
