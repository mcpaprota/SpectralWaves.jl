
"""
    update_bbc_sle!(A′, A″, w′, Ψ̂′, Ψ̂″, Ψ̃′, Ψ̃″, β̂, κ, κ′, ℐ, F, M, msg_flag)

Compute coefficients `A′`, `A″`, and `w′` for the bottom boundary condition
system of linear equations (SLE).

Modified in-place variables:
- `A′` are constant coefficients of the system of linear equations and
- `A″` are coefficients of the system of linear equations,
- `w′` are coefficients corresponding to the linear wavemaker term.

"""
function update_bbc_sle!(A′, A″, w′, Ψ̂′, Ψ̂″, Ψ̃′, Ψ̃″, β̂, κ, κ′, ℐ, F, M, msg_flag)
    N, _ = convolution_range(0, M, ℐ)
    b̃ = complex(zeros(N))
    b̂ = complex(zeros(N))
    A′[:,:] = zeros(2ℐ + 1, 2ℐ + 1)
    A″[:,:] = diagm(ones(2ℐ + 1))
    β̃ = im * κ .* β̂
    w′[:] = (2im * κ′ * β̃ + β̂ ^ 2)
    _, r1 = convolution_range(1, M, ℐ)
    msg_flag ? progr = Progress(M, barglyphs=BarGlyphs("[=> ]"), color=:black) : nothing
    for m in 0:M-1
        _, r = convolution_range(m, M, ℐ)
        b̃[r] = β̂ ^ m * β̃ / F[m+1]
        B̃ = toeplitz(b̃[r1])
        b̂[r] = β̂ ^ m * β̂ / F[m+2]
        B̂ = toeplitz(b̂[r1])
        A′[:,:] += B̃ .* transpose(Ψ̃′[:, m+1]) - B̂ .* transpose(Ψ̂′[:, m+2]) # SLE constant coefficient matrix
        A″[:,:] += - B̃ .* transpose(Ψ̃″[:, m+1]) + B̂ .* transpose(Ψ̂″[:, m+2]) # SLE coefficient matrix
        msg_flag && next!(progr)
    end
    msg_flag && finish!(progr)
    return nothing
end

"""
    nonlinear_kfsbc_correction(η̂, ϕ̂, ψ̂, Φ̂′, Φ̂″, Φ̃′, Φ̃″, κ, κ′, ℐ, F, M)

Calculate nonlinear correction `δη̇` to kinematic free-surface boundary condition.

"""
function nonlinear_kfsbc_correction(η̂, ϕ̂, ψ̂, Φ̂′, Φ̂″, Φ̃′, Φ̃″, κ, κ′, ℐ, F, M, ξ, ℓ)
    N, r0 = convolution_range(0, M, ℐ)
    δη̇ = complex(zeros(N))
    η̃ = im * κ .* η̂
    Φ̃ = Φ̃′ .* ϕ̂ + Φ̃″ .* ψ̂
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
    nonlinear_dfsbc_correction(η̂, ϕ̂, ϕ̇, ψ̂, ψ̇, Φ̇′, Φ̇″, Φ̂′, Φ̂″, Φ̃′, Φ̃″, κ′, ℐ, F, M, ξ, ζ, ℓ, d)

Calculate nonlinear correction `δϕ̇` to dynamic free-surface boundary condition.

"""
function nonlinear_dfsbc_correction(η̂, ϕ̂, ϕ̇, ψ̂, ψ̇, Φ̇′, Φ̇″, Φ̂′, Φ̂″, Φ̃′, Φ̃″, κ′, ℐ, F, M, ξ, ζ, ℓ, d)
    N, r0 = convolution_range(0, M, ℐ)
    δϕ̇ = complex(zeros(N))
    Φ̇ = Φ̇′ .* ϕ̇ + Φ̇″ .* ψ̇
    Φ̃ = Φ̃′ .* ϕ̂ + Φ̃″ .* ψ̂
    Φ̂ = Φ̂′ .* ϕ̂ + Φ̂″ .* ψ̂
    Φ̇[ℐ+1, 2] += 2ζ * d / ℓ
    Φ̃[:, 1] -= 2im * ξ / ℓ * κ′
    Φ̂[ℐ+1, 2] += 2ξ / ℓ
    if M > 1
        Φ̇[ℐ+1, 3] += 2ζ / ℓ
    end
    for m in 0:M-1
        Φ² = complex(zeros(4ℐ + 1))
        for n in 0:m
            Φ² += binomial(m, n) * (Φ̃[:, n+1] * Φ̃[:, m-n+1] + Φ̂[:, n+1] * Φ̂[:, m-n+1])
        end
        _, r = convolution_range(m + 1, M, ℐ)
        δϕ̇[r] += η̂ ^ m / F[m+1] * (η̂ * Φ̇[:, m+2] / (m + 1) + Φ² / 2)
    end
    return δϕ̇[r0]
end

"""
    time_integration_coeffs(O::Integer)

Calculate Adams-Bashforth-Moulton time-stepping scheme coefficients for a given order `O`.

Output is a tuple of two vectors with Adams-Bashforth and Adams-Moulton coefficients,
respectively.

"""
function time_integration_coeffs(O::Integer)
    O == 1 && return [1.0], [1.0]
    O == 2 && return [3, -1] / 2, [1, 1] / 2
    O == 3 && return [23, -16, 5] / 12, [5, 8, -1] / 12
    return [55, -59, 37, -9] / 24, [9, 19, -5, 1] / 24
end


"""
    solve_problem!(p::Problem; msg_flag = true)

Calculate solution coefficients `η̂`, `η̇`, `ϕ̂`, `ϕ̇`, `ψ̂`, `ψ̇` of the wave problem.

Modified in-place variables:
- `η̂` are free-surface elevation potential amplitudes (m),
- `η̇` are free-surface vertical velocity amplitudes (m/s),
- `ϕ̂` are flat-bottom velocity potential amplitudes (m²/s),
- `ϕ̇` are flat-bottom acceleration potential amplitudes (m²/s²),
- `ψ̂` are uneven-bottom velocity potential amplitudes (m²/s),
- `ψ̇` are uneven-bottom acceleration potential amplitudes (m²/s²),

Input variables:
- `β̂` are bottom-surface elevation amplitudes (m),
- `β̇` are bottom-surface vertical velocity amplitudes (m/s),
- `p̂` are surface pressure head amplitudes (m),
- `κ` are wave numbers (rad/m),
- `𝒯` are hyperbolic tangent lookup values,
- `𝒮` are hyperbolic secant lookup values,
- `ℐ` is the number of harmonics,
- `M_s` is the order of nonlinear free-surface boundary condition,
- `M_b` is the order of nonlinear bottom boundary condition,
- `Δt` is the time step (s),
- `O` is the order of the time-stepping scheme,
- `N` is the number of time steps,
- `χ` is wavemaker paddle displacement (m),
- `ξ` is wavemaker paddle velocity (m/s),
- `ζ` is wavemaker paddle acceleration (m/s²),
- `ℓ` is the fluid domain length (m),
- `d` is the water depth (m).

Keyword arguments:
- `msg_flag` controls whether progress messages are printed.

"""
function solve_problem!(p::Problem; msg_flag = true)
    ℓ, d = p.ℓ, p.d
    ℐ, Δt, N, O, M_s, M_b, F = p.ℐ, p.Δt, p.N, p.O, p.M_s, p.M_b, p.F
    ϕ̂, ϕ̇, ψ̂, ψ̇ = p.ϕ̂, p.ϕ̇, p.ψ̂, p.ψ̇
    η̂, η̇, p̂ = p.η̂, p.η̇, p.p̂
    β̂, β̇ = p.β̂, p.β̇
    ξ, ζ = p.ξ, p.ζ
    κ, 𝒯, 𝒮 = p.κ, p.𝒯, p.𝒮
    static_bottom = p.static_bottom
    # initialize auxiliary variables
    c_ab, c_am = time_integration_coeffs(O)
    κ′ = @.  1 / κ * (κ ≠ 0)
    κ″ = @.  1 / κ^2 * (κ ≠ 0)
    # initialize nonlinear bottom boundary condition if necessary
    if M_b > 0
        msg_flag && println("Initializing system of equations...")
        A′, A″, Ψ̂′, Ψ̂″, Ψ̃′, Ψ̃″, w′ = init_nonlinear_bottom_boundary_condition(κ, 𝒯, 𝒮, ℐ, M_b)
        if static_bottom
            msg_flag && println("Factorization of the system for static bottom...")
            update_bbc_sle!(A′, A″, w′, Ψ̂′, Ψ̂″, Ψ̃′, Ψ̃″, β̂[:, 1], κ, κ′, ℐ, F, M_b, msg_flag)
            A″ = factorize(A″)
        end
    end
    # initialize nonlinear free-surface boundary conditions if necessary
    if M_s > 0
        Φ̇′, Φ̇″, Φ̂′, Φ̂″, Φ̃′, Φ̃″ = init_nonlinear_surface_boundary_condition(κ, 𝒯, 𝒮, ℐ, M_s)
        # define short calls to nonlinear correction functions
        δη̇(n) = nonlinear_kfsbc_correction(η̂[:, n], ϕ̂[:, n], ψ̂[:, n], Φ̂′, Φ̂″, Φ̃′, Φ̃″, κ, κ′, ℐ, F, M_s, ξ[n], ℓ)
        δϕ̇(n) = nonlinear_dfsbc_correction(η̂[:, n], ϕ̂[:, n], ϕ̇[:, n], ψ̂[:, n], ψ̇[:, n], Φ̇′, Φ̇″, Φ̂′, Φ̂″, Φ̃′, Φ̃″, κ′, ℐ, F, M_s, ξ[n], ζ[n], ℓ, d)
    end
    msg_flag && println("Time marching...")
    msg_flag ? progress_bar = Progress(N, barglyphs=BarGlyphs("[=> ]"), color=:black, showspeed=true) : nothing
    # start time-marching loop
    for n in O:N+O-1
        j = 0
        # initial first guess of acceleration potential amplitudes
        if (n > O) && (M_s > 0)
            @views ϕ̇[:, n] = ϕ̇[:, n-1]
            @views ψ̇[:, n] = ψ̇[:, n-1]
        end
        # initialize loop for iterative solver to wave problem
        while j < J
            # apply dynamic free-surface boundary condition
            if M_s == 0
                @views ϕ̇[:, n] = -g * (η̂[:, n] + p̂[:, n]) + 2ζ[n] * κ″ / ℓ
            else
                @views ϕ̇[:, n] = -g * (η̂[:, n] + p̂[:, n]) + 2ζ[n] * κ″ / ℓ - δϕ̇(n)
            end
            ϕ̇[ℐ+1, n] = -g * η̂[ℐ+1, n] - ζ[n] * (d^2 / ℓ - ℓ / 12)
            # apply Adams-Bashforth predictor
            @views ϕ̂[:, n+1] = ϕ̂[:, n] + Δt * sum(c_ab[i] * ϕ̇[:, n+1-i] for i in 1:O)
            # apply nonlinear bottom boundary condition if necessary
            if M_b > 0
                if !static_bottom
                    update_bbc_sle!(A′, A″, w′, Ψ̂′, Ψ̂″, Ψ̃′, Ψ̃″, β̂[:, n+1], κ, κ′, ℐ, F, M_b, false)
                end
                b = A′ * ϕ̂[:, n+1] + β̇[:, n+1] - ξ[n+1] / ℓ * w′[ℐ+1:3ℐ+1]
                ψ̂[:, n+1] = A″ \ b
            end
            # apply kinematic free-surface boundary condition
            @views η̇[:, n+1] = @. κ * 𝒯 * ϕ̂[:, n+1] + ψ̂[:, n+1] * 𝒮
            η̇[ℐ+1, n+1] += 2ξ[n+1] * d / ℓ
            if M_s > 0
                @views η̇[:, n+1] -= δη̇(n+1)
            end
            η̂ₚ = η̂[:, n+1]
            # apply Adams-Moulton corrector
            @views η̂[:, n+1] = η̂[:, n] + Δt * sum(c_am[i] * η̇[:, n+2-i] for i in 1:O)
            if M_s > 0
                # check accuracy of the solution
                general_error(η̂ₚ, η̂[:, n+1]) < ϵ ? break : j += 1
                # apply central difference scheme
                n == O ? ψ̇[:, n] = ψ̂[:, n+1] / Δt : ψ̇[:, n] = (ψ̂[:, n+1] - ψ̂[:, n-1]) / 2Δt
                # check simulation blow-up
                if !isfinite(norm(η̂[:, n+1]))
                    msg_flag && finish!(progress_bar)
                    msg_flag && println("Finished unsuccessfully!")
                    return false
                end
            else
                break
            end
        end
        msg_flag && next!(progress_bar)
    end
    msg_flag && finish!(progress_bar)
    msg_flag && println("Finished successfully!")
    return true
end
