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

# infix version of convolve
Base.:*(a::Vector{<:Number}, b::Vector{<:Number}) = convolve(a, b)

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

Compute range `r` of indices of a central part of convolution vector
corresponding to the `m`-th order of convolution,
where `M` is the maximum order of convolution,
`N` is the length of convolution vector corresponding to the `M`-th order of convolution,
and `ℐ` is the number of harmonics.

# Examples
```julia-repl
julia> m = 1; M = 3; ℐ = 5; c_range = convolution_range(m, M, ℐ)
11:31
```
"""
function convolution_range(m::Integer, M::Integer, ℐ::Integer)
    N = (M + 1) * 2ℐ + 1
    r = (N+1)÷2-(m+1)*ℐ:(N+1)÷2+(m+1)*ℐ
    return N, r
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

"""
    general_error(a::Vector{<:Number}, b::Vector{<:Number})

Compute a general error between vectors `a` and `b` as a `relative_error(a, b)`
or an `absolute_error(a, b)` for `norm(a)` not equal or equal to 0, respectively.

"""
function general_error(a::Vector{<:Number}, b::Vector{<:Number})
    norm(a) == 0 && return absolute_error(a, b)
    return relative_error(a, b)
end


"""
    factorial_lookup(n::Integer)

Compute lookup table for factorials of numbers from 0 to `n`.

# Examples
```julia-repl
julia> factorial_lookup(3)
4-element Vector{Float64}:
  1.0
  1.0
  2.0
  6.0
```
"""
function factorial_lookup(n::Integer)
    F = convert(Vector{Float64}, factorial.(big.(0:n)))
    return F
end

"""
    inverse_fourier_transform(f̂::Vector{<:Number}, ω::AbstractRange{<:Number}, x::Number)

Compute `f(x)` using expansion amplitudes `f̂` and eigenvalues `ω`.

# Examples
```julia-repl
julia> f̂ = [0.5, 0, 0.5]; ω = -1:1; x = π; inverse_fourier_transform(f̂, ω, x)
-1.0
```
"""
function inverse_fourier_transform(f̂::Vector{<:Number}, ω::AbstractRange{<:Number}, x::Number)
    f = real(sum(f̂ .* exp.(im * ω * x)))
    return f
end

"""
    fourier_transform(f::Vector{<:Number}, ω::Number, x::AbstractRange{<:Number})

Compute `f̂(ω)` using function values `f` at points `x`.

# Examples
```julia-repl
julia> f = [-1, 0, 1, 0, -1]; ω = 1; x = range(-π, π, length = 5);
julia> fourier_transform(f, ω, x)
0.5 + 0.0im
```
"""
function fourier_transform(f::Vector{<:Real}, ω::Number, x::AbstractRange{<:Number})
    k = ones(Integer, length(f))
    k[[1, end]] = [2, 2]
    f̂ = sum(f ./ k .* exp.(-im * ω * x)) * (x[2] - x[1]) / (x[end] - x[1])
    return f̂
end
