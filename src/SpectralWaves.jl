module SpectralWaves

using LinearAlgebra, ProgressMeter

include("algebra.jl")
include("constants.jl")
include("problem.jl")
include("solvers.jl")
include("waves.jl")
include("bottoms.jl")
include("output.jl")

# inner functions exported only for testing for now
export Problem
export g, ϵ, J
export convolve, toeplitz, convolution_power, relative_error, absolute_error
export convolution_range, factorial_lookup, inverse_fourier_transform, fourier_transform
export init_problem, solve_problem!
export linear_regular_wave!, linear_wavemaker!
export water_surface, bottom_surface
export surface_bump!, bottom_bump!

end
