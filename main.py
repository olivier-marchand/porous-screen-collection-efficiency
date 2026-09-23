#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
##################################################
##################################################
# Porous screen collection efficiency
### ---------------------------- ###
# Version 1.0.0
### ---------------------------- ###
# Supplementary material for Inertial impact of particles on porous screen, 
# available at <https://doi.org/10.1103/ch74-gbk3>.
# Co-authors include Olivier C. Marchand, Christophe Josserand and
# Camille Duprat.
#
# For description, we refer to the README.md file.
#
# Contact: olivier.marchand@univ-grenoble-alpes.fr
### ---------------------------- ###
# MIT License
#
# Copyright (c) 2026 Olivier C. Marchand
##################################################
##################################################
"""

__authors__ = 'Olivier C. Marchand, Christophe Josserand, Camille Duprat'
__contact__ = 'Olivier C. Marchand <olivier.marchand@univ-grenoble-alpes.fr>'
__copyright__ = 'Copyright (c) 2026 Olivier Claude Marchand'
__license__ = 'MIT Licence'
__date__ = '21/09/2026'
__version__ = '1.0.0'

#################################################
############### Module importation ##############
#################################################

"""
# Import aerodynamic_efficiency.py and capture_efficiency.py, 
# must be in the same directory as this file.
"""

from aerodynamic_efficiency import *
from capture_efficiency import *


#################################################
############### Capture efficiency ##############
#################################################

# --------- # Dimensionless parameters

rp_adim = 0.1
St_d = 0.1

# --------- # Function

eta_c = efficiency_eta_c(rp_adim,St_d)
print('# --------- #')
print('Capture efficiency: ',eta_c)

#%%

#################################################
############ Aerodynamic efficiency #############
#################################################

# --------- # Dimensionless parameters

s = 0.5
Re_d = 100
St_L = 0.1

# --------- # Function

eta_a = efficiency_eta_a(s,Re_d,St_L)
print('# --------- #')
print('Aerodynamic efficiency: ',eta_a)

#%% Alternative use

#################################################
#### Computation from dimensional parameters ####
#################################################

# --------- # Dimensional parameters

v0 = 3 # (m/s)
d = 0.001 # (m)
L = 1 # (m)
nu_f = 1.5*10**(-5) # (m^2/s)
mu_f = 1.8*10**(-5) # (Pa.s)
rho_p = 1000 # (kg/m^3)
d_p = 10**(-5) # (m)
S = 0.3 # (m^2)

# --------- # Computation of the velocity close to the screen

s = S/L**2
print("s = ",s)

Re_d = v0*d/nu_f
print("Re_d = ",Re_d)

v_n = (1-1/2*x_3D_final(s, Re_d))*v0

# --------- # Cpmputation of the other dimensionless parameters

rp_adim = d_p/d
print("rp_adim = ",rp_adim)

St_d = 1/18*rho_p*d_p**2*v_n/(mu_f*d)
print("St_d = ",St_d)

St_L = 1/18*rho_p*d_p**2*v0/(mu_f*L)
print("St_L = ",St_d)

eta_c = efficiency_eta_c(rp_adim,St_d)
print('# --------- #')
print('Capture efficiency: ',eta_c)

eta_a = efficiency_eta_a(s,Re_d,St_L)
print('# --------- #')
print('Aerodynamic efficiency: ',eta_a)

# --- #
