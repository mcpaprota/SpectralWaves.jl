module SpectralWaves

using LinearAlgebra, ProgressMeter

include("algebra_functions.jl")
include("constants.jl")
include("problem.jl")
include("solvers.jl")
include("waves.jl")

# inner functions exported only for testing for now
export g, ϵ, J
export convolve, toeplitz, convolution_power, relative_error, absolute_error
export convolution_range, factorial_lookup, inverse_fourier_transform, fourier_transform
export init_problem, solve_problem!
export linear_regular_wave!, linear_wavemaker!

end
