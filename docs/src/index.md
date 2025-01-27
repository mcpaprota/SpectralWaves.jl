```@meta
CurrentModule = SpectralWaves
```

# SpectralWaves.jl

A Fourier Galerkin method solution to nonlinear waves propagating over topography.

## Overview

[`SpectralWaves.jl`](https://github.com/mcpaprota/SpectralWaves.jl) is a Julia package for simulation of nonlinear waves propagating over arbitrary topography under potential flow assumptions. The solution is derived using a Fourier Galerkin spectral method in terms of amplitudes of free-surface elevation and velocity potential, while inverse Fourier transform is used to get a phase-resolved wave field. Four wave generation mechanisms are supported - initial conditions, linear wavemaker forcing, pressure forcing (to be implemented), moving bottom.

## Wave problem

We consider waves propagating over arbitrary bottom topography in a periodic fluid domain of length $\ell$ and characteristic depth $d$ (corresponding to still water level). A Cartesian coordinate system is used to define fluid elements along horizontal $x$-axis coinciding with undisturbed free surface and upward-pointing and vertical $z$-axis. Undulating free surface is described by means of $\eta(x, t)$, while bottom topography is considered as fluctuations $\beta(x, t)$ around $-d$. The general scheme is presented in _figure 1_ below.

```@example
include("../../examples/docs_figure_1.jl") # hide
fig # hide
```
_Figure 1: Propagation of waves over topography in a periodic fluid domain._

According to potential flow assumptions (irrotational flow of an inviscid and incompressible fluid), we define velocity vector field $\mathbf{v}(x, z, t) = \nabla\varPhi(x, z, t)$ and formulate a our boundary-value problem in a following way (_table 1_)

_Table 1: Initial boundary-value problem of waves propagating over topography._

| Equation | Region | Description |
|:--------|:------|:-----------|
| $$ \frac{\partial^2\varPhi}{\partial x^2} + \frac{\partial^2\varPhi}{\partial z^2} = 0 $$ | $$ -d + \beta \leq z \leq \eta $$ | Laplace's equation |
| $$ \frac{\partial\eta}{\partial t} + \frac{\partial\eta}{\partial x}\frac{\partial\varPhi}{\partial x} - \frac{\partial\varPhi}{\partial z} = 0 $$ | $$ z = \eta $$ | Kinematic free-surface boundary condition |
| $$ \frac{\partial\varPhi}{\partial t} + g\eta + \frac{1}{2}\left(u^2 + w^2 \right) = 0 $$ | $$ z = \eta $$ | Dynamic free-surface boundary condition |
| $$ \frac{\partial\beta}{\partial t} + \frac{\partial\beta}{\partial x}\frac{\partial\varPhi}{\partial x} - \frac{\partial\varPhi}{\partial z} = 0 $$ | $$ z = \beta - d $$ | Kinematic bottom boundary condition |
| $$ \varPhi(x + \ell, z, t) = \varPhi(x, z, t) $$ | $$ x = 0, \ell $$ | Periodic lateral boundary condition |

where $u = \partial\varPhi / \partial x$ and $w  = \partial\varPhi / \partial z$ are horizontal and vertical velocity components, respectively, while a gravitational acceleration $g\approx 9.81\,\mathrm{m/s}^2$.

## Spectral solution

We use spectral expansions of $\varPhi$, $\eta$, and $\beta$, while additionally $\varPhi$ is decomposed into parts: $\phi$ - satisfying homogeneous problem of waves propagating over horizontal bottom and $\psi$ - satisfying a corrugated bottom correction. The total velocity potential $\varPhi$ satisfies Laplace equation. In _Table 2_, we provide spectral expansion formulas of $\varPhi$, $\phi$, $\psi$, $\eta$, and $\beta$ along with velocity components $u$ and $w$.

_Table 2: Spectral expansion formulas._

| Equation | Description |
|--------|----|
| $$ \varPhi(x, z, t) = \phi(x, z, t) + \psi(x, z, t) $$ | Total velocity potential |
| $$ \phi(x, z, t) = \sum_{i=-I}^{I}\hat{\phi}_i(t)\frac{\cosh\kappa_i(z+d)}{\cosh\kappa_id}\mathrm{e}^{\mathrm{i}\kappa_ix} $$ | Flat-bottom velocity potential |
| $$ \psi(x,z,t) = \hat{\psi_0}(t)z + \sum_{i=-I\,\wedge\,i\ne0}^{I}\hat{\psi}_i(t)\frac{\sinh\kappa_iz}{\kappa_i\cosh\kappa_id}\mathrm{e}^{\mathrm{i}\kappa_ix}$$ | Corrugated bottom velocity potential |
| $$ \eta(x, t) = \sum_{i=-I}^{I}\hat{\eta}_i(t)\mathrm{e}^{\mathrm{i}\kappa_ix} $$ | Free-surface elevation |
| $$ \beta(x, t) = \sum_{i=-I}^{I}\hat{\beta}_i(t)\mathrm{e}^{\mathrm{i}\kappa_ix} $$ | Bottom topography |
| $$ u(x, z, t) = \sum_{i=-I\,\wedge\,i\ne0}^{I}\mathrm{i}\frac{\hat{\phi}_i(t)\kappa_i\cosh\kappa_i(z+d) + \hat{\psi}_i(t)\sinh\kappa_iz}{\cosh\kappa_id}\mathrm{e}^{\mathrm{i}\kappa_ix}$$ | Horizontal velocity component|
| $$ w(x, z, t) = \hat{\psi_0}(t)+ \sum_{i=-I\,\wedge\,i\ne0}^{I}\frac{\hat{\phi}_i(t)\kappa_i\sinh\kappa_i(z+d) + \hat{\psi}_i(t)\cosh\kappa_iz}{\cosh\kappa_id}\mathrm{e}^{\mathrm{i}\kappa_ix}$$ | Vertical velocity component|