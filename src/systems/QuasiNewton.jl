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

# Example
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
# Example
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
# Example
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
# Example
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
    fdJ(X, FX = F(X), h = 10^-4, n = size(X,1))

(Forward) Finite Difference approximation of the Jacobian.
Requires starting vector X. Optional h and n.
FX is optional in case F(X) was already computed for efficiency.
# Example
```julia-repl
julia> fdJ([1.5;2.0])
2×2 Array{Float64,2}:
 3.0001  -1.0
 1.0     -4.0001
```
"""
function fdJ(X,FX = F(X),h::Float64 = 10^-4,n = size(X,1))
    H = h * I(n)
    J = X #placeholder to set dims
    for i = 1:n
        J = hcat(J, F(X + H[:,i]) - FX)
    end
    J = h^-1 * J[1:end, 2:end] #removes placeholder and divide by h
end

"""
    FiniteDiff(X, h = 10^-4, eps = 10^-6)

Generalized Finite-Difference Method for a system of nonlinear equations in [`F`](@ref).
Requires starting vector X. Optional h and eps.
# Example
```julia-repl
julia> FiniteDiff([1.5;2.0])
2-element Array{Float64,1}:
 1.6180340699246842
 1.6180341416981485
```
"""
function FiniteDiff(X,h::Float64 = 10^-4, eps::Float64 = 10^-6)
    n = size(X,1)
    #println(X) #uncomment to see itr
    FX = F(X)
    while(norm(FX) > eps)
        J = fdJ(X,FX,h,n)
        X = X - J \ FX
        #println(X) #uncomment to see itr
        FX = F(X)
    end
    X
end

"""
    Broyden(X, eps = 10^-6)

"Good" Broyden Method for a system of nonlinear equations in [`F`](@ref).
Requires starting vector X. Optional eps. Initial J approximated by differencing.
# Example
```julia-repl
julia> Broyden([1.5;2.0])
2-element Array{Float64,1}:
 1.6180340721138766
 1.6180340507162783
```
"""
function Broyden(X, eps::Float64 = 10^-6)
    #println(X) #uncomment to see itr
    FX = F(X)
    invJ = inv(fdJ(X,FX))
    while(norm(FX) > eps)
        Xold = X
        X = X - invJ*FX
        #println(X) #uncomment to see itr
        deltaX = X - Xold
        FXold = FX
        FX = F(X)
        deltaF = FX - FXold
        trans = transpose(deltaX)
        invJ = invJ + (deltaX - invJ*deltaF)/(trans*invJ*deltaF)*trans*invJ
    end
    X
end

"""
Main
"""

println(Secant(transpose([1.0 1.0;1.0 2.0;1.5 2.0])))
println(Steffensen([1.5;2.0]))
println(FiniteDiff([1.5;2.0]))
println(Broyden([1.5;2.0]))
