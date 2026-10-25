# Constants for the SpectralWaves.jl package

"""
    g

Gravitational acceleration in metres per second squared.
"""
const g::Float64 = 9.81 # gravitational acceleration (m/s^2)

"""
    ϵ

Relative error tolerance used for convergence of the free-surface elevation.
"""
const ϵ::Float64 = 1e-12 # relative error tolerance for free-surface elevation

"""
    J

Maximum number of iterations used by the nonlinear solver.
"""
const J::Integer = 10000 # maximum number of iterations
