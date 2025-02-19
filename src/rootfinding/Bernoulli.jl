using Polynomials # just used for printing

"""
# Bernoulli's Method and Quotient-Difference Method

This file contains examples of using both methods for finding roots.

# Main Functions
- `aitken`: Aitken's Δ² Method.
- `bernolli`: Finds the dominant root.
- `bernoulli_small`: Finds the smallest root.
- `qd`: Finds all roots (assuming real and distrinct).

# Dependencies
- `Polynomial.jl` package used for printing polynomials

# Author
Oscar Veliz

# Date
February 14, 2025 
"""

"""
    aitken(a, b, c)

Aitken's Δ² Method to accelerate convergence of a sequence.

# Arguments
- `a`: The first term in the sequence.
- `b`: The second term in the sequence.
- `c`: The third term in the sequence.

# Returns
Predicted later value in sequence.

# Example
```julia
julia> aitken(1.0, 2.0, 1.5)
1.625
```
"""
function aitken(a, b, c)
    return a - (b - a)^2 / (c - 2 * b + a)
end


"""
    bernoulli(a::Vector{<:Number}, with_a::Bool=false, eps::Float64=1e-8, max_iter::Int=60) -> Union{Float64, NaN}

Perform Bernoulli's Method to find the dominant root of a polynomial.

# Arguments
- `a::Vector{<:Number}`: A vector of numbers representing the coefficients of a polynomial. Example `[-1.0,-1.0,1.0]` is x²-x-1.
- `with_aitken::Bool`: A boolean flag to determine whether to use Aitken's Method. Default is `false`.
- `eps::Float64`: A small value used as the convergence threshold. Default is `1e-8`.
- `max_iter::Int`: The maximum number of iterations allowed. Default is `60`.

# Returns
- `Union{Float64, NaN}`: The root or `NaN` if not found.

# See also
- `aitken`
- `bernoulli_small`

# Example
```julia
julia> bernoulli([-1, -1, 1])
1.6180339901755971
```
"""
function bernoulli(a::Vector{<:Number}, with_aitken=false, eps=1e-8, max_iter=60)
    # println(Polynomial(a)) # uncomment to see polynomial
    n = length(a)
    x = zeros(max_iter)
    q = zeros(max_iter)
    x[n-1] = 1.0
    t = n
    j = 1
    while t < max_iter
        if with_aitken && t != n && mod(j, 4) == 0 && j > 3
            q[j] = aitken(q[j-1], q[j-2], q[j-3])
            x[t] = q[j] * x[t-1]
        else # not aitken
            x[t] = -sum(a[k] * x[t-mod(n - k, n)] for k in 1:n-1) / a[n]
            q[j] = x[t] / x[t-1]
        end

        if j > 2 * n && abs(q[j] - q[j-1]) <= eps # 2 cycles and small enough steps
            # println(j, " iterations") # uncomment to see iterations
            # println("q ",q[1:j]) # uncomment to see q
            # println("x ",x[1:t]) # uncomment to see x            
            return q[j]
        end
        j += 1
        t += 1
    end
    # println(j-1, " iterations") # uncomment to see iterations
    # println("x ", x) # uncomment to see x
    # println("q ", q[1:j-1]) # uncomment to see q
    return NaN # in case a result not found
end


"""
    bernoulli_small(a::Vector{<:Number}, with_a::Bool=false, eps::Float64=1e-8, max_iter::Int=60) -> Union{Float64, NaN}

Perform Bernoulli's Method to find the smallest root of a polynomial. Invokes `bernoulli` on an inverted polynomial.

# Arguments
- `a::Vector{<:Number}`: A vector of numbers representing the coefficients of a polynomial. Example `[-1.0,-1.0,1.0]` is x²-x-1.
- `with_aitken::Bool`: A boolean flag to determine whether to use Aitken's Method. Default is `false`.
- `eps::Float64`: A small value used as the convergence threshold. Default is `1e-8`.
- `max_iter::Int`: The maximum number of iterations allowed. Default is `60`.

# Returns
- `Union{Float64, NaN}`: The root or `NaN` if not found.

# See also
- `aitken`
- `bernoulli`

# Example
```julia
julia> bernoulli_small([-1, -1, 1])
-0.6180339882053251
```
"""
function bernoulli_small(a::Vector{<:Number}, with_aitken=false, eps=1e-8, max_iter=60)
    c = reverse(a)
    c ./= c[end] # normalize - divide every element by leading coefficient
    # println(c) # uncomment to see new coefficients
    return 1.0 / bernoulli(c, with_aitken, eps, max_iter)
end


"""
    qd(a::Vector{<:Number}, eps::Float64=1e-8, max_iter::Int=60) -> Vector{Float64}

Perform Quotient-Difference Algorithm to find all the root of a polynomial. This is the unstable version that only finds real distinct roots.

# Arguments
- `a::Vector{<:Number}`: A vector of numbers representing the coefficients of a polynomial. Example `[-1.0,-1.0,1.0]` is x²-x-1.
- `eps::Float64`: A small value used as the convergence threshold. Default is `1e-8`.
- `max_iter::Int`: The maximum number of iterations allowed. Default is `60`.

# Returns
- `Vector{Float64}`: The roots found.

# See also
- `bernoulli`

# Example
```julia
julia> bernoulli_small([-1, -1, 1])
2-element Vector{Float64}:
  1.618033985017358
 -0.618033993331295
```
"""
function qd(a::Vector{<:Number}, eps=1e-8, max_iter=60)
    # println(Polynomial(a)) # uncomment to see polynomial
    n = length(a)
    x = zeros(max_iter)
    q = zeros((max_iter, n - 1))
    e = zeros((max_iter, n - 1))
    x[n-1] = 1.0
    t = n
    j = 1
    # initial q -- Bernoulli's Method
    while t < max_iter
        x[t] = -sum(a[k] * x[t-mod(n - k, n)] for k in 1:n-1) / a[n]
        q[j, 1] = x[t] / x[t-1]
        if j > n && abs(q[j, 1] - q[j-1, 1]) <= eps # 2 cycles and small enough steps
            # println(j, " iterations") # uncomment to see iterations
            break
        end
        j += 1
        t += 1
    end
    max_iter = j - 1

    if (n == 2) # single root don't try QD
        return [q[max_iter, 1]]
    end

    # QD Table
    for k in 2:n-1
        for i in 2:max_iter
            e[i-1, k] = q[i, k-1] - q[i-1, k-1] + e[i, k-1]
        end
        for i in 2:max_iter
            q[i, k] = e[i+1, k] / e[i, k] * q[i+1, k-1]
        end
    end
    # println(j-1, " iterations") # uncomment to see iterations
    # println("x ", x[1:max_iter]) # uncomment to see x
    # println("e = ", e[1:max_iter,1:n-1]) # uncomment to see e
    # println("q = ", q[1:max_iter,1:n-1]) # uncomment to see q

    # Check for instability
    roots = zeros(n - 1)
    roots[1] = q[max_iter, 1]
    for i in 2:n-1
        last_step = 10000
        for k in 2:max_iter-1
            next_step = abs(q[k, i] - q[k+1, i])
            if next_step > last_step
                roots[i] = q[k, i]
                break
            end
            last_step = next_step
        end
    end
    return roots
end


"""
    main()

Examples of using Bernoulli's Method and QD to find roots. Invoked when program is executed.
"""
function main()
    # coeff = [3, 1.0]
    coeff = [-1.0, -1.0, 1.0]
    # coeff = [-2310, 727, 382, -72, -8, 1.0]
    println("Largest root of ", Polynomial(coeff), " ≈ ", bernoulli(coeff, false, 1e-6))
    println("Smallest root of ", Polynomial(coeff), " ≈ ", bernoulli_small(coeff, false, 1e-6))
    println("All roots of ", Polynomial(coeff), " ≈ ", qd(coeff))

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
