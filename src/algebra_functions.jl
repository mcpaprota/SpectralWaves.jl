"""
    convolve(a::Vector{T}, b::Vector{T}) where T

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

Base.:^(a::Vector{<:Number}, n::Integer) = convolution_power(a, n) # infix version of convolution n-th power

"""
    toeplitz(a::Vector{T}) where T

Transform vector `a` to a Toeplitz matrix.

# Examples
```julia-repl
julia> a = [1, 2, 3]; toeplitz(a)
3×3 Matrix{Int64}:
 2  1
 3  2
```
"""
function toeplitz(a::Vector{T}) where T<:Number
    n = Int((length(a) - 1) ÷ 2)
    B = zeros(T, n+1, n+1)
    for i in 0:n, j in 0:n
        B[i+1, j+1] = a[i-j+n+1]
    end
    return B
end
