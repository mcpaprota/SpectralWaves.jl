module SpectralWaves

using LinearAlgebra, ProgressMeter

include("algebra_functions.jl")
include("constants.jl")

# inner functions exported only for testing for now
export convolve, toeplitz, convolution_power, relative_error, absolute_error

end
