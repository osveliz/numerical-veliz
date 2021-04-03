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
x - y^2 + 1

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
2-element Array{Float64,1}:
 1.6180339736070495
 1.6180339507485015
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
    steffJ(X, FX = F(X), n = size(X,1))

Traub-Steffensen approximation of the Jacobian.
Requires starting vector X.
FX is optional in case F(X) was already computed for efficiency.
# Examples
```julia-repl
julia> steffJ([1.5;2.0])
2×2 Array{Float64,2}:
 2.25  -1.0
 1.0   -2.5
```
"""
function steffJ(X,FX = F(X),n = size(X,1))
    H = Diagonal(FX)
    J = X #placeholder to set dims
    for i = 1:n
        J = hcat(J, F(X + H[:,i]) - FX)
    end
    J = J[1:end, 2:end] * inv(H) #removes placeholder before mult
end

"""
    Steffensen(X, eps=10^-6)

Traub-Steffensen Method for a system of nonlinear equations in [`F`](@ref).
Requires starting vector X.
# Examples
```julia-repl
julia> Steffensen([1.5;2.0])
2-element Array{Float64,1}:
 1.6180339736216265
 1.6180339387805316
```
"""
function Steffensen(X, eps::Float64 = 10^-6)
    n = size(X,1)
    #println(X) #uncomment to see itr
    FX = F(X)
    while(norm(FX) > eps)
        J = steffJ(X,FX,n)
        X = X - J \ FX
        #println(X) #uncomment to see itr
        FX = F(X)
    end
    X
end


"""
Coming Soon ... additional Quasi-Newton Methods
"""

"""
Main
"""

println(Secant(transpose([1.0 1.0;1.0 2.0;1.5 2.0])))
println(Steffensen([1.5;2.0]))
