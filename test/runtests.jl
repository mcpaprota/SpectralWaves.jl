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
    @test SpectralWaves.fourier_transform([0, 1, 0], 0, 0:2) == 0.5 + 0.0im
end
