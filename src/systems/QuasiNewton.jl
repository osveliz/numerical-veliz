using LinearAlgebra
"""
Quasi-Newton Methods for Solving Systems of Nonlinear Equations
@author Oscar Veliz
"""

"""
    F(X)

System of nonlinear equations on vector X ([x;y]).
Overwrite for your own testing, currently set to:
x^2 - y - 1
y - x^2 + 1

# Examples
```julia-repl
julia> F([1.0; 2.0])
[-2.0, -2.0]
```
"""
F(X) = [X[1]^2 - X[2] - 1; X[1] - X[2]^2 + 1]


"""
    Secant(X, eps=10^-6)

Wolfe-Bittner Secant Method for the system of nonlinear equations in [`F`](@ref).
Requires `d+1` starting vectors in X.
# Examples
```julia-repl
julia> Secant(transpose([1.0 1.0;1.0 2.0;1.5 2.0]))
[1.6180339736070495, 1.6180339507485015]
```
"""
function Secant(X, eps::Float64 = 10^-6)
    n = size(X,1)
    FN = mapslices(F, X; dims=1)
    FX = FN[:,end]
    #println(X) #uncomment to see starting points
    while norm(FX) > eps
        A = [ones(1,n+1); FN]
        b = zeros(n+1,1)
        b[1,1] = 1
        p = A \ b
        Xbar = X * p
        #print(Xbar) #uncomment to see each itr
        X = hcat(X, Xbar)
        X = X[1:end,2:end]
        FX = F(Xbar)
        FN = hcat(FN,FX)
        FN = FN[1:end,2:end]
    end
    X[:,end]
end


"""
Coming Soon ... additional Quasi-Newton Methods
"""

"""
Main
"""

println(Secant(transpose([1.0 1.0;1.0 2.0;1.5 2.0])))
