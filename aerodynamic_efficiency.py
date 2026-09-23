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
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
import time
from scipy.interpolate import splprep, splev

#################################################
############ Potential flow variables ###########
#################################################

def pressure_drop(s_var,Re):
    """ 
    Pressure drop accross a screen composed of fibres, 
    Re is the local Reynolds number with the 
    velocity v_n at the vicinity of the screen.
    
    Note: this function must be adapted to the screen type. 
    """
    return((11/Re+0.8/np.log(Re+1.25)+0.055*np.log(Re))*(1-1/(1-s_var)**2))

def func_eqn(x,Red_var,s_var):
    """ 
    Equation corresponding to the conservation of the momentum 
    in a control volume around the screen.
    """
    gamma0 = 0.0998 # Geometric parameter, see Marchand et al. (2024).
    k = pressure_drop(s_var,Red_var*(1-1/2*x[0])) 
    return [-1/8*x[0]**4*k+x[0]**2*(8*gamma0 + k - 2) - 4*x[0] - 2*k]

def x_3D_final(s,Red):
    """ 
    Compute the value of the dimensionless potential source omega.
    """
    # ---------- # Different starting estimate depending on Red.
    if s<=0.2:
        if Red<10:
            x0=1.99
        else:
            x0 = 0.3
    if s>0.2 and s<0.4:
        if Red<10:
            x0=1.99
        else:
            x0 = 1.0
    if s>=0.4 and s<0.6:
        if Red<10:
            x0=1.99
        else:
            x0 = 1.3
    if s>=0.6 and s<0.8:
        if Red<10:
            x0=1.999
        else:
            x0 = 1.8
    if s>=0.8 and s<0.95:
        x0 = 1.95
    if s>=0.95:
        x0 = 1.999
    if Red<10:
        if s>0.9:
            x0 = 1.999
    # ---------- # Resolution of the equation.
    root = fsolve(func_eqn, [x0], args=(Red,s))
    x = root[0] # The first root is the correct one.
    return(x)

#################################################
################ Velocity field #################
#################################################

def FI(x,y,z,u,v):
    return(1/(4*np.pi)*np.log(np.sqrt(z**2+(x-u)**2+(v-y)**2)+v-y))

def FJ(x,y,z,u,v):
    return(1/(4*np.pi)*np.log(np.sqrt(z**2+(x-u)**2+(v-y)**2)+u-x))

def FK(x,y,z,u,v):
    return(1/(4*np.pi)*np.arctan((x-u)*(v-y)/(-z*np.sqrt(z**2+(x-u)**2+(v-y)**2))))

def vIx(x,y,z,omega):
    return(omega*(FI(x,y,z,-1/2,-1/2)-FI(x,y,z,1/2,-1/2)-FI(x,y,z,-1/2,1/2)+FI(x,y,z,1/2,1/2)))

def vIy(x,y,z,omega):
    return(omega*(FJ(x,y,z,-1/2,-1/2)-FJ(x,y,z,1/2,-1/2)-FJ(x,y,z,-1/2,1/2)+FJ(x,y,z,1/2,1/2)))

def vIz(x,y,z,omega):
    return(1+omega*(FK(x,y,z,-1/2,-1/2)-FK(x,y,z,1/2,-1/2)-FK(x,y,z,-1/2,1/2)+FK(x,y,z,1/2,1/2)))

#################################################
##### Resolution of the equation of motion ######
#################################################

def func(t,y,St,omega):
    return([y[3],y[4],y[5],-1/St*(y[3]-vIx(y[0],y[1],y[2],omega)),-1/St*(y[4]-vIy(y[0],y[1],y[2],omega)),-1/St*(y[5]-vIz(y[0],y[1],y[2],omega))])

def sol_p(St,tfadim,y0,omega):
    """ 
    Solve the motion equation of a particle 
    in a potential flow around a 3D square porous screen.
    """
    t0adim = 0 # s
    T = np.linspace(t0adim,tfadim,10000)
    sol = solve_ivp(func, [t0adim,tfadim], y0, t_eval = T, method='Radau', args=(St,omega))
    return(sol)

#################################################
####### Search for limiting trajectories ########
#################################################

def func_impact(omega,x0adim,y0adim,z0adim,St_Lval,rp_Ladim):
    """ 
    Check if there is an impact for the given initial positions.
    """
    y0init = [x0adim,y0adim,z0adim,vIx(x0adim,y0adim,z0adim,omega),vIy(x0adim,y0adim,z0adim,omega),vIz(x0adim,y0adim,z0adim,omega)]
    Sol = sol_p(St_Lval,63,y0init,omega)
    impact_val = 0
    Z = Sol.y[2]
    vxval = None
    vyval = None
    vzval = None
    for i in range(0,len(Z)):
        if Z[i] >= -rp_Ladim and abs(Sol.y[0][i])<=1/2 and abs(Sol.y[1][i])<=1/2:
            impact_val = 1
            vxval = Sol.y[3][i]
            vyval = Sol.y[4][i]
            vzval = Sol.y[5][i]
            #print("x = ",Sol.y[0][i])
            #print("y = ",Sol.y[1][i])
            break
        else:
            None
    return([impact_val,vxval, vyval, vzval])


def func_curve_traj_tang(omega,St_Lval,rp_Ladim,z0adim,nr):
    """
    Compute the upstream positions of the particules 
    following the tangent trajectories.
    """
    # ---------- # Precision paremeters.
    r0prec = 0
    epsilon0 = 10**(-3)
    # ---------- # Initial position of the particules.
    Theta = [i*np.pi/(4*(nr-1)) for i in range(0,nr)]
    r0_init = 1/2*np.sqrt(2)
    R = []
    Vx_impadim = []
    Vy_impadim = []
    Vz_impadim = []
    for i in range(0,nr):
        r0 = r0_init
        r0prec = 0
        theta0 = Theta[i]
        #print("################")
        #print("New trajectory!")
        #print("theta = ",theta0*360/(2*np.pi))
        # ---------- # Initialization
        epsilon = 1
        detect_0 = 0
        impact = 0
        nb_iteration = 0
        # ---------- # Iteration until the precision is meet.
        while epsilon >= epsilon0:
            nb_iteration = nb_iteration + 1
            #print("---------------")
            #print("Iteration : ",nb_iteration)
            #print("r0 = ",r0)
            x0adim = r0*np.sin(theta0)
            y0adim = r0*np.cos(theta0)
            Fimp = func_impact(omega,x0adim,y0adim,z0adim,St_Lval,rp_Ladim)
            impact = Fimp[0]
            epsilon =  abs((r0 - r0prec)/2)
            #print("epsilon = ",epsilon)
            rstack = r0
            if impact == 0: 
                r0 = r0 - epsilon
                #print("No impact!")
            if impact == 1:
                detect_0 = 1
                r0 = r0 + epsilon
                #print("Impact!")
                vxval = Fimp[1]
                vyval = Fimp[2]
                vzval = Fimp[3]
            r0prec = rstack
        # ---------- #
        if detect_0==0:
            # No impact.
            vxval = None
            vyval = None
            vzval = None
        # ---------- # The impact velocity is also returned.
        R.append(r0)
        Vx_impadim.append(vxval)
        Vy_impadim.append(vyval)
        Vz_impadim.append(vzval)
    return([R,Theta,Vx_impadim,Vy_impadim,Vz_impadim])

#################################################
### Computation of the aerodynamic efficiency ###
#################################################

def generation_data(svar,Red_var,St_Lval,rp_Ladim,z0adim):
    """
    Intermediate function to compute 
    the aerodynamic efficiency.
    """
    nr = 10 # Number of discretized points (over 1/8 of the total screan surface).
    # -------------------------- #
    xvar = x_3D_final(svar, Red_var)
    #print("vn_adim = ",1-1/2*xvar)
    omega = xvar
    E = (1-omega/2)/(1+omega/2)
    
    #print("omega = ", omega)
    #print("E = ", E)
    
    Result = func_curve_traj_tang(omega,St_Lval,rp_Ladim,z0adim,nr)
    return(Result)

def complete_curve(Theta0,R0,nr):
    Theta_total = []
    R_reversed = []
    R_total = []
    k = nr-1
    for i in range(0,nr):
        Theta_total.append(Theta0[i])
        R_total.append(R0[i])
        R_reversed.append(R0[i])
    R_reversed.reverse()
    for i in range(1,nr):
        k = k + 1
        Theta_total.append(k*np.pi/4*1/(nr-1))
        R_total.append(R_reversed[i])
    for i in range(1,nr):
        k = k + 1
        Theta_total.append(k*np.pi/4*1/(nr-1))
        R_total.append(R0[i])
    for i in range(1,nr):
        k = k + 1
        Theta_total.append(k*np.pi/4*1/(nr-1))
        R_total.append(R_reversed[i])
    for i in range(1,nr):
        k = k + 1
        Theta_total.append(k*np.pi/4*1/(nr-1))
        R_total.append(R0[i])
    for i in range(1,nr):
        k = k + 1
        Theta_total.append(k*np.pi/4*1/(nr-1))
        R_total.append(R_reversed[i])
    for i in range(1,nr):
        k = k + 1
        Theta_total.append(k*np.pi/4*1/(nr-1))
        R_total.append(R0[i])
    for i in range(1,nr-1):
        k = k + 1
        Theta_total.append(k*np.pi/4*1/(nr-1))
        R_total.append(R_reversed[0:-1][i])
    return([R_total,Theta_total])

def area_polygon(x,y):
    '''
    Area of a polygon defined 
    by len(x) points with coordinates (x,y),
    using the cross product formula.
    '''
    return(1/2*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1))))

def efficiency_eta_a(svar,Red_var,St_Lval):
    """
    Compute the aerodynamic efficiency.
    ----------
    Parameters
    ----------
    svar : float
        Screen solidity between 0 and 1.
    Red_var : float
        Local Reynolds number, typically between 1 and 1000.
    St_Lval : float
        Global Stokes number.
    ----------
    Note: rp_adim is set to zero by default only for the computation 
    of the aerodynamic efficiency.
    """
    start = time.time()
    # ---------- #
    X1 = []
    Y1 = []
    # ---------- #
    rp_adim = 0
    z0adim = -60 # Initial particle released at z=60L upstream, z0adim = z0/L.
    # ---------- #
    Data = generation_data(svar,Red_var,St_Lval,rp_adim,z0adim)
    Result = [Data[0],Data[1]]
    # ---------- # Reconstructs the total interception surface contour
    COMP = complete_curve(Result[1],Result[0],len(Result[0]))
    for i in range(0,len(COMP[0])):
        X1.append(COMP[0][i]*np.sin(COMP[1][i]))
        Y1.append(COMP[0][i]*np.cos(COMP[1][i]))
    # ---------- # Reconstructs the interception surface area from the discretized contour
    tck, u = splprep([X1, Y1], s=0, per=1)
    new_points = splev(u, tck)
    eta_a = area_polygon(new_points[0],new_points[1])
    # ---------- #
    end = time.time()
    #print("Computational time: ", end - start)
    return(eta_a)

