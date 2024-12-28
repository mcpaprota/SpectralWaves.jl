"""
    convolve(a, b) or a *b

Compute direct convolution of two vectors `a` and `b`.

# Examples
```julia-repl
julia> a = [1, 2]; b = [1, 1]; convolve(a,b)
3-element Vector{Int64}:
 1
 3
 2
```
"""
function convolve(a, b)
    n = length(a)
    m = length(b)
    c = zero([a; b[1:end-1]])
    for i in 1:n
        for j in 1:m
            c[i+j-1] += a[i] * b[j]
        end
    end
    return c
end

Base.:*(a, b) = convolve(a, b) # infix version of convolution

"""
    convolution_power(a, n) or a ^ n

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
function convolution_power(a, n::Integer)
    m = length(a)
    if n == 0
        return a[1] * 0 + 1
    elseif n == 1
        return a
    else
        a * convolution_power(a, n - 1)
    end
end

Base.:^(a, n) = convolution_power(a, n) # infix version of convolution n-th power

"""
    toeplitz(a)

Transform vector `a` to Toeplitz matrix.

# Examples
```julia-repl
julia> a = [1, 2, 3, 4, 5]; toeplitz(a)
3×3 Matrix{Int64}:
 3  2  1
 4  3  2
 5  4  3
```
"""
function toeplitz(a)
    n = Int((length(a) - 1) ÷ 2)
    B = similar(a, n + 1, n + 1) * 0
    for i in 0:n, j in 0:n
        B[i+1, j+1] = a[i-j+n+1]
    end
    return B
end
