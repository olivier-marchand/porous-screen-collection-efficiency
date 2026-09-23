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

import numpy as np
#import time
from scipy.integrate import solve_ivp

#################################################
##### Resolution of the equation of motion ######
#################################################

def func(t,y,St):
    return([y[1],y[3]**2*y[0]-1/St*(y[1] - (1 - 1/y[0]**2)*np.cos(y[2])),y[3],-2*y[1]*y[3]/y[0]-1/St*(y[3]+1/y[0]*(1+1/y[0]**2)*np.sin(y[2]))])

def sol_p(St,tfadim,y0):
    """ 
    Solve the equation of motion of a particle in a potential flow around a cylinder.
    """
    t0adim = 0
    T = np.linspace(t0adim,tfadim,10000)
    sol = solve_ivp(func, [t0adim,tfadim], y0, t_eval = T, method='Radau', args=(St,))
    return(sol)

#################################################
######## Search for limiting trajectory #########
#################################################

def func_impact(z0,x0,St_val,rp_adim):
    """ 
    Check if there is an impact for the given initial positions.
    """
    vz0_val = 1 
    vx0_val = 0
    tfadim = 22
    r01 = np.sqrt(z0**2+x0**2)
    theta01 = np.arctan(x0/z0) + np.pi
    y01 = [r01,vz0_val*np.cos(theta01)+vx0_val*np.sin(theta01),theta01,1/r01*(vx0_val*np.cos(theta01)-vz0_val*np.sin(theta01))]
    S1 = sol_p(St_val,tfadim,y01)
    # ---------- # Condition d'impact.
    impact_val = 0
    for i in range(0,len(S1.y[0])):
        if S1.y[0][i] > 1 + rp_adim:
            None
        else:
            impact_val = 1
            break
    return(impact_val)

def func_x0_traj_tang(St_val,rp_adim,z0):
    """ 
    Compute the tangent trajectory by dichotomy method.
    """
    x0_init = 1
    x0 = x0_init
    x0prec = 0
    xstack = x0_init
    epsilon = 1
    epsilon0 = 10**(-4)
    impact = 0
    nb_iteration = 0
    while epsilon >= epsilon0:
        nb_iteration = nb_iteration + 1
        #print("---------------")
        #print("Iteration: ",nb_iteration)
        impact = func_impact(z0,x0,St_val,rp_adim)
        epsilon =  abs((x0 - x0prec)/2)
        #print(epsilon)
        xstack = x0
        if impact == 0: 
            x0 = x0 - epsilon
            #print("No impact!")
        if impact == 1:
            x0 = x0 + epsilon
            #print("Impact!")
        x0prec = xstack
    return([x0,x0prec])

#################################################
### Computation of the aerodynamic efficiency ###
#################################################

def find_xc(St_val,rp_adim,z0):
    """ 
    Compute the ordinate at -infty for the tangent trajectory 
    """
    x0 = func_x0_traj_tang(St_val,rp_adim,z0)[0]
    xc = x0*(1-1/(x0**2+z0**2))
    #print("###########")
    #print("St = ",St_val)
    #print("xc = ",xc)
    #print("###########")
    return(xc)

def efficiency_eta_c(rp_adim,Std):
    #start = time.time()
    """
    Compute the capture efficiency.
    ----------
    Parameters
    ----------
    rp_adim : float
        Particle to fibre diameter ratio.
    Std : float
        Local Stokes number.
    ----------
    Note: For very high Stokes number, z0adim should be increased.
    Possibly, the value tfadim in func_impact should be adapted.
    """
    z0adim = -20 # Initial particle released at z=20d upstream, z0adim = 2*z0/d.
    eta_c = find_xc(Std,rp_adim,z0adim)
    #end = time.time()
    #print("Computational time: ", end - start)
    return(eta_c)