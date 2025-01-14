
"""
    update_bbc_sle!(A′, A″, β̃, β̂, Ψ̂′, Ψ̂″, Ψ̃′, Ψ̃″, κ, ℐ, F, M)

Compute coefficients `A′`, `A″`, and `β̃` for the bottom boundary condition
system of linear equations.

"""
function update_bbc_sle!(A′, A″, β̃, β̂, Ψ̂′, Ψ̂″, Ψ̃′, Ψ̃″, κ, ℐ, F, M)
    N, _ = convolution_range(0, M, ℐ)
    a = complex(zeros(N))
    A′[:] = zeros(2ℐ + 1, 2ℐ + 1)
    A″[:] = diagm(ones(2ℐ + 1))
    β̃[:] = im * κ .* β̂
    _, r1 = convolution_range(1, M, ℐ)
    for m in 0:M-1
        _, r = convolution_range(m, M, ℐ)
        a[r] = β̂ ^ m * β̃ / F[m+1]
        B̃ = toeplitz(a[r1])
        a[r] = β̂ ^ m * β̂ / F[m+1]
        B̂ = toeplitz(a[r1])
        A′[:] = @. A′ + B̃ * transpose(Ψ̃′[:, m+1]) - B̂ * transpose(Ψ̂′[:, m+2]) # SLE constant coefficient matrix
        A″[:] = @. A″ - B̃ * transpose(Ψ̃″[:, m+1]) + B̂ * transpose(Ψ̂″[:, m+2]) # SLE coefficient matrix
    end
end

"""
    nonlinear_kfsbc_correction(η̂, ϕ̂, ψ̂, Φ̂′, Φ̂″, Φ̃′, Φ̃″, κ, κ′, ℐ, F, M)

Calculate nonlinear correction `δη̇` to kinematic free-surface boundary condition.

"""
function δη̇(η̂, ϕ̂, ψ̂, Φ̂′, Φ̂″, Φ̃′, Φ̃″, κ, κ′, ℐ, F, M, ξ, ℓ)
    N, r0 = convolution_range(0, M, ℐ)
    δη̇ = complex(zeros(N))
    η̃ = im * κ .* η̂
    Φ̃ = Φ̃′ .* ϕ̃ + Φ̃″ .* ψ̃
    Φ̂ = Φ̂′ .* ϕ̂ + Φ̂″ .* ψ̂
    Φ̃[:, 1] -= 2im * ξ / ℓ * κ′
    Φ̂[ℐ+1, 2] += 2ξ / ℓ
    for m in 0:M-1
        _, r = convolution_range(m + 1, M, ℐ)
        δη̇[r] += η̂ ^ m / F[m+1] * (η̃ * Φ̃[:, m+1] - η̂ * Φ̂[:, m+2] / (m + 1))
    end
    return δη̇[r0]
end

"""
    nonlinear_dfsbc_correction(η̂, ϕ̂, ψ̂, Φ̂′, Φ̂″, Φ̃′, Φ̃″, κ′, ℐ, F, M, ξ, ζ, ℓ, d)

Calculate nonlinear correction `δϕ̇` to dynamic free-surface boundary condition.

"""
function δϕ̇(η̂, ϕ̂, ψ̂, Φ̂′, Φ̂″, Φ̃′, Φ̃″, κ′, ℐ, F, M, ξ, ζ, ℓ, d)
    N, r0 = convolution_range(0, M, ℐ)
    δϕ̇ = complex(zeros(N))
    Φ̇ = Φ̇′ .* ϕ̃ + Φ̇″ .* ψ̃
    Φ̃ = Φ̃′ .* ϕ̃ + Φ̃″ .* ψ̃
    Φ̂ = Φ̂′ .* ϕ̂ + Φ̂″ .* ψ̂
    Φ̇[ℐ+1, 2] += 2ζ * d / ℓ
    Φ̃[:, 1] -= 2im * ξ / ℓ * κ′
    Φ̂[ℐ+1, 2] += 2ξ / ℓ
    if M > 1
        Φ̇[ℐ+1, 3] += 2ζ / ℓ
    end
    for m in 0:M-1
        Φ² = complex(zeros(4ℐ + 1))
        for n in 0:M
            Φ² = binomial(m, n) * (Φ̃[:, n+1] * Φ̃[:, m-n+1] + Φ̂[:, n+1] * Φ̂[:, m-n+1])
        end
        _, r = convolution_range(m + 1, M, ℐ)
        δϕ̇[r] += η̂ ^ m / F[m+1] * (η̂ * Φ̇[:, m+2] / (m + 1) + Φ²)
    end
    return δϕ̇[r0]
end
