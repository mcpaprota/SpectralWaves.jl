"""
    convolve(a::Vector{Ta}, b::Vector{Tb}) where {Ta<:Number, Tb<:Number}

Compute direct convolution of vectors `a` and `b`.

# Examples
```julia-repl
julia> a = [1, 2]; b = [1, 1]; convolve(a, b)
3-element Vector{Int64}:
 1
 3
 2
```
"""
function convolve(a::Vector{Ta}, b::Vector{Tb}) where {Ta<:Number, Tb<:Number}
    n = length(a)
    m = length(b)
    T = promote_type(Ta, Tb)
    c = zeros(T, n + m - 1)
    for i in 1:n
        for j in 1:m
            c[i+j-1] += a[i] * b[j]
        end
    end
    return c
end

Base.:*(a::Vector{<:Number}, b::Vector{<:Number}) = convolve(a, b) # infix version of convolve

"""
    convolution_power(a::Vector{<:Number}, n::Integer)

Compute convolution `n`-th power of vector `a`.

# Examples
```julia-repl
julia> a = [1, 2, 3]; convolution_power(a, 3)
7-element Vector{Int64}:
  1
  6
 21
 44
 63
 54
 27
```
"""
function convolution_power(a::Vector{<:Number}, n::Integer)
    if n == 0
        return a[1] * 0 + 1
    elseif n == 1
        return a
    else
        a * convolution_power(a, n - 1)
    end
end

# infix version of convolution n-th power
Base.:^(a::Vector{<:Number}, n::Integer) = convolution_power(a, n)

"""
    convolution_range(m::Integral, M::Integral, n::Integer)

Compute range of indices of a central part of convolution vector corresponding to
'm'-th order of convolution, where `M` is the maximum order of convolution and `n` is
the number of harmonics.

# Examples
```julia-repl
julia> m = 1; M = 3; n = 5; c_range = convolution_range(m, M, n)
11:31
```
"""
function convolution_range(m::Integer, M::Integer, n::Integer)
    N = (M + 1) * 2n + 1
    r = (N+1)÷2-(m+1)*n:(N+1)÷2+(m+1)*n
    return r
end

"""
    toeplitz(a::Vector{T}) where T

Transform vector `a` to a Toeplitz matrix.

# Examples
```julia-repl
julia> a = [1, 2, 3]; toeplitz(a)
2×2 Matrix{Int64}:
 2  1
 3  2
```
"""
function toeplitz(a::Vector{T}) where T<:Number
    n = (length(a) - 1) ÷ 2
    B = zeros(T, n+1, n+1)
    for i in 0:n, j in 0:n
        B[i+1, j+1] = a[i-j+n+1]
    end
    return B
end

"""
    relative_error(a::Vector{<:Number}, b::Vector{<:Number})

Compute relative error between vectors `a` and `b`.

# Examples
```julia-repl
julia> a = [2, 2]; b = [1, 1]; relative_error(a, b)
0.5
```
"""
function relative_error(a::Vector{<:Number}, b::Vector{<:Number})
    δ = norm(a - b) / norm(a)
    return δ
end

"""
    absolute_error(a::Vector{<:Number}, b::Vector{<:Number})

Compute absolute error between vectors `a` and `b`.

# Examples
```julia-repl
julia> a = [1, 2]; b = [1, 1]; absolute_error(a, b)
1.0
```
"""
function absolute_error(a::Vector{<:Number}, b::Vector{<:Number})
    ϵ = norm(a - b)
    return ϵ
end
