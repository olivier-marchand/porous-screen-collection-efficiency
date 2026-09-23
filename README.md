# Porous screen collection efficiency

This folder contains codes which supplement the article *Inertial impact of particles on porous screen* available at <https://doi.org/10.1103/ch74-gbk3>.

## Description

The folder contains three files that enable the computation of the capture efficiency and aerodynamic efficiency associated to the collection of particles by a porous square screen composed of fibres using the following dimensional parameters:

- `v0` (m/s): the usptream fluid velocity.
- `d` (m): the fibre diameter.
- `L` (m): the length of the square screen.
- `nu_f` (m^2/s): the fluid kinematic viscosity.
- `mu_f` (Pa.s): the fluid dynamic viscosity.
- `rho_p` (kg/m^3): the density of the particle.
- `d_p` (m): the diameter of the particle.
- `S` (m^2): the area of the solid part of the screen.

Alternatively, the user can set the following input dimensionless parameters:

- `s`: the solidity of the screen computed from `S/L**2`.
- `Re_d`: the local Reynolds number associated to the pore scale, computed from `Re_d=v0*d/nu_f`.
- `St_d`: the local Stokes number computed from `1/18*rho_p*d_p**2*v_n/(mu_f*d)` with `v_n` the velocity at the vicinity of the screen (provided by the model if not imposed by the user).
- `St_L`: the global Stokes number computed from `1/18*rho_p*d_p**2*v0/(mu_f*L)`.
- `rp_adim`: the ratio between the particle diameter and the fibre diameter, thus computed from `r_p_adim=d_p/d`.

The present computation of the trajectory of the particle contains several assumptions that are discussed in the mentioned paper. In particular the local Reynolds number `Re_d` should not be much lower than 1. The computation of the capture efficiency assumes that the particle is captured once it comes into contact with a solid part of the screen. Therefore, this part of the code should be used carefully for non-liquid particles.

## Getting Started

### Dependencies

The last version of the code requires

- Python 3.11.14

with installed packages

- numpy 2.4.2
- scipy 1.16.3

### Executing Program

An example of executing code is provided in the file `main.py`. The files `capture_efficiency.py` and `aerodynamic_efficiency.py` must be placed in the same folder than the executing program if `main.py` is used (otherwise the path must be adapted to import the modules).

## Authors

Contributors include Olivier C. Marchand, Christophe Josserand and Camille Duprat.

Questions can be adressed to Olivier C. Marchand at <olivier.marchand@univ-grenoble-alpes.fr>.

## Version History

- 1.0.0
    - Initial release.

## Licence

This project is licensed under the MIT Licence - see the LICENCE.md file for details.
