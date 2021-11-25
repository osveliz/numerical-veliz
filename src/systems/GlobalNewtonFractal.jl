using LinearAlgebra
using Plots
gr()

"""
Global Newton Fractal for Solving Systems of Nonlinear Equations
@author Oscar Veliz
"""

"""
    F(X)

System of nonlinear equations on vector X ([x;y]).
Overwrite for your own testing, currently set to:
x^2 - y - 1
x - y^2 + 1

For a system version of z^3-1 use F(X) = [X[1]^3-3*X[1]*X[2]^2-1;3*X[1]^2*X[2]-X[2]^3]

# Example
```julia-repl
julia> F([1.0; 2.0])
2-element Array{Float64,1}:
 -2.0
 -2.0
```
"""
F(X) = [X[1]^2-X[2]-1;X[1]-X[2]^2+1]

"""
    J(X)

Jacobian matrix for F given vector X ([x;y]).
Overwrite for your own testing, currently set to:
[2x -1;1 -2y]

For a system version of z^3-1 use J(X) = [3*X[1]^2-3*X[2]^2 -6*X[1]*X[2];6*X[1]*X[2] 3*X[1]^2-3*X[2]^2]

# Example
```julia-repl
julia> J([1.0; 2.0])
2×2 Array{Float64,2}:
 2.0  -1.0
 1.0  -4.0
```
"""
J(X) = [2*X[1] -1;1 -2*X[2]]

"""
    Newton(X, eps=10^-6)

Newton's Method for the system of nonlinear equations in [`F`](@ref). Returns the number of iterations.
# Example
```julia-repl
julia> Newton([1.0;2.0])
5
```
"""
function Newton(X, eps::Float64=10^-6)
    i = 1
    FX = F(X)
    while norm(FX) > eps && i < 50
        X = X - J(X) \ FX
        FX = F(X)
        i = i + 1
    end
    i
end

"""
    GlobalNewton(X, eps=10^-6)

Global Newton's Method for the system of nonlinear equations in [`F`](@ref). Returns the number of iterations.
# Example
```julia-repl
julia> GlobalNewton([1.0;2.0])

```
"""
function GlobalNewton(X, eps::Float64=10^-6)
    FX = F(X)
    nF = norm(FX)
    i = 1
    JX = J(X)
    delta = 1000
    while norm(delta) > eps && nF > eps
        a = 1.0
        delta = JX \ FX
        T = X - a * delta
        FT = F(T)
        nFT = norm(FT)
        while nFT > nF
            a = a * 0.5
            T = X - a * delta
            FT = F(T)
            nFT = norm(FT)
        end
        delta = X - T
        X = T
        FX = FT
        nF = nFT
        JX = J(X)
        i = i + 1
    end
    i
end

"""
    Fractal(useGlobal=true, eps=10^-3)

Function for creating fractals. Call with false for normal newton fractal. Uses [`F`](@ref) and [`J`](@ref).
Saves the file to `globalnewton.svg` or `newton.svg` depending on `useGlobal`.
# Example
```julia-repl
julia> Fractal()
```
"""
function Fractal(useGlobal::Bool=true, step::Float64=.001)
    domain = -3.6:step:3.6
    range = -2.0:step:2.0
    data = rand(length(range), length(domain))
    i = 1
    for x in domain
        j = 1
        for y in range
            X = [x;y]
            try
                data[j,i] = if useGlobal GlobalNewton(X) else Newton(X) end
            catch
                data[j,i] = 50
            end
            j = j + 1
        end
        if i > length(domain)
            i = 1
        else
            i = i + 1
        end
    end
    heatmap(domain,
    range, data,
    c = cgrad([:black, :red, :yellow, :white]),
    xlabel = "x", ylabel = "y",
    title = if useGlobal "Global Newton Fractal" else "Newton Fractal" end,
    fmt = :svg)
    savefig(if useGlobal "globalnewton.svg" else "newton.svg" end)
end

"""
Main
"""

Fractal()
