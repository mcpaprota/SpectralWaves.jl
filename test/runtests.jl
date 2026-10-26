using SpectralWaves
using Test

@testset "SpectralWaves.jl" begin
    # convolve function test
    @test SpectralWaves.convolve([1, 2], [1, 1]) == [1, 3, 2]
    @test SpectralWaves.convolve([1, 2], [1, 1im]) == [1 + 0im, 2 + 1im, 0 + 2im]
    @test [1, 2] * [1, 1] == [1, 3, 2]
    # convolution_power function test
    @test SpectralWaves.convolution_power([1, 2], 0) == 1
    @test SpectralWaves.convolution_power([1, 2], 1) == [1, 2]
    @test SpectralWaves.convolution_power([1, 2, 3], 2) == [1, 4, 10, 12, 9]
    @test [1, 2] ^ 0 == 1
    @test [1, 2] ^ 1 == [1, 2]
    # toeplitz function test
    @test SpectralWaves.toeplitz([1, 2, 3]) == [2 1; 3 2]
    # relative_error function test
    @test SpectralWaves.relative_error([2, 2], [1, 1]) == 0.5
    # absolute_error function test
    @test SpectralWaves.absolute_error([1, 2], [1, 1]) == 1.0
    # convolution_range function test
    @test SpectralWaves.convolution_range(1, 3, 5) == (41, 11:31)
    # factorial_lookup function test
    @test SpectralWaves.factorial_lookup(3) == [1.0, 1.0, 2.0, 6.0]
    # inverse_fourier_transform function test
    @test SpectralWaves.inverse_fourier_transform([0.5, 0, 0.5], -1:1, π) == -1
    # fourier_transform function test
    @test SpectralWaves.fourier_transform([-1, 0, 1, 0, -1], 1, range(-π, π, length = 5)) == 0.5 + 0.0im
    # solving Problem test floats
    ℓ = 1.0
    d = 1.0
    ℐ = 1
    t = range(0, 1, length = 10)
    O = 4
    M_s = 0
    M_b = 0
    static_bottom = true
    problem = Problem(ℓ, d, ℐ, t; O = O, M_s = M_s, M_b = M_b, static_bottom = static_bottom)
    @test solve_problem!(problem; msg_flag = false)
    # solving Problem test integers and range of floats
    ℓ = 1
    d = 1
    ℐ = 1
    t = range(0, 1, length = 10)
    O = 4
    M_s = 0
    M_b = 0
    static_bottom = true
    problem = Problem(ℓ, d, ℐ, t; O = O, M_s = M_s, M_b = M_b, static_bottom = static_bottom)
    @test solve_problem!(problem; msg_flag = false)
    # solving Problem test integers
    ℓ = 1
    d = 1
    ℐ = 1
    t = 1:10
    O = 4
    M_s = 0
    M_b = 0
    static_bottom = true
    problem = Problem(ℓ, d, ℐ, t; O = O, M_s = M_s, M_b = M_b, static_bottom = static_bottom)
    @test solve_problem!(problem; msg_flag = false)

    # convolve_dsp agrees with convolve
    a, b = [1.0, 2.0, 3.0], [0.5, -1.0]
    @test convolve_dsp(a, b) ≈ convolve(a, b)
    # general_error switches between absolute and relative error
    @test SpectralWaves.general_error([0, 0], [3, 4]) == absolute_error([0, 0], [3, 4])
    @test SpectralWaves.general_error([2, 2], [1, 1]) == relative_error([2, 2], [1, 1])
    # fourier_transform followed by inverse_fourier_transform recovers a cosine
    x = range(-π, π, length = 201)
    f̂ = [fourier_transform(cos.(x), k, x) for k in -2:2]
    @test inverse_fourier_transform(f̂, -2:2, 0.7) ≈ cos(0.7) atol = 1e-3
end

@testset "Problem construction" begin
    p = Problem(2π, 1.0, 3, range(0, 1, length = 5); O = 3, M_s = 2, M_b = 1, static_bottom = false)
    @test p.N == 4
    @test p.Δt == 0.25
    @test size(p.η̂) == (7, 7)
    @test size(p.β̂) == (7, 7)
    @test p.κ == -3:3
    @test p.F == [1.0, 1.0, 2.0]
    @test size(Problem(1, 1, 2, 1:5).β̂) == (5, 1)
end

@testset "Bottom shapes" begin
    p = Problem(4.0, 1.0, 8, range(0, 1, length = 5))
    bottom_bump!(p, 0.1, 1.0, 0.5)
    @test bottom_surface(p, 0.5) > bottom_surface(p, 1.5)
    bottom_step!(p, 0.2)
    @test bottom_surface(p, 1.9) ≈ 0.2 atol = 0.05
    @test bottom_surface(p, 0.0) ≈ 0.0 atol = 0.05
    @test bottom_surface(p, -1.0) ≈ bottom_surface(p, 1.0)
    bottom_slope!(p, 0.2)
    @test bottom_surface(p, 1.9) ≈ 0.2 atol = 0.05
    @test bottom_surface(p, 0.0) ≈ 0.0 atol = 0.05
    x = range(-2, 2, length = 101)
    bottom_vector!(p, x, fill(0.3, length(x)))
    @test bottom_surface(p, 0.0) ≈ 0.3 atol = 1e-6
    pm = Problem(4.0, 1.0, 4, range(0, 1, length = 5); static_bottom = false)
    moving_bottom_bump!(pm, 0.1, 1.0, 0.5)
    @test any(!iszero, pm.β̂)
    @test any(!iszero, pm.β̇)
    @test bottom_surface(pm, 0.0, pm.O) isa Real
end

@testset "Waves and outputs" begin
    ℓ, d, ω = 2π, 1.0, 2.0
    p = Problem(ℓ, d, 4, range(0, 1, length = 5))
    linear_regular_wave!(p, 0.1, ω)
    @test water_surface(p, 0.0, 1) ≈ 0.05 * 2 / 2 atol = 1e-12
    @test water_surface(p, π, 1) ≈ -water_surface(p, 0.0, 1)
    @test water_velocity(p, 0.0, 0.0, 1, :x) isa Real
    @test water_velocity(p, 0.0, 0.0, 1, :z) isa Real
    @test water_velocity(p, 0.0, 0.0, 1, :y) === nothing

    pb = Problem(ℓ, d, 4, range(0, 1, length = 5))
    surface_bump!(pb, 0.1, 1.0, 1.0)
    @test water_surface(pb, 1.0, 1) > water_surface(pb, 4.0, 1)

    pp = Problem(ℓ, d, 4, range(0, 1, length = 5); static_bottom = false)
    moving_pressure_bump!(pp, 0.1, 1.0, 0.5)
    @test surface_pressure(pp, 0.0, pp.O) > 0

    pw = Problem(ℓ, d, 4, range(0, 2, length = 21))
    linear_wavemaker!(pw, 0.1, 1.0, 2.0, 1)
    @test pw.χ[pw.O] == 0
    @test any(!iszero, pw.ξ) && any(!iszero, pw.ζ)
end

@testset "Linear solve with nonlinear orders" begin
    p = Problem(2π, 1.0, 4, range(0, 0.5, length = 6); M_s = 1, M_b = 1)
    linear_regular_wave!(p, 0.01, 1.0)
    bottom_step!(p, 0.1)
    @test solve_problem!(p; msg_flag = false)
    @test all(isfinite, abs.(p.η̂))
end
