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
end
