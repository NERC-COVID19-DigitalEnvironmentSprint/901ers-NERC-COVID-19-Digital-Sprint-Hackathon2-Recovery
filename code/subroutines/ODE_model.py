# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 10:13:19 2020

@author: Paul
"""

import numpy as np
import matplotlib.pyplot as plt


def ODE_model(t,dt,u,pars,Y0,Lockdown=False):
    
    """
    Solving single time step of Covid model using Forward Euler
    
    Inputs:
        
        - t, the current time (units:yr)
        
        - dt, the current time-step (units:yr)
        
        - Array of current variables:
            - Y, the current GDP (units: £billiion)
            - eta, the current carbon intensity
            - E, the current emissions
            - R, the current reproduction value
            - m, the current mobility
            - mu, the current decarbonisation
            
        - Array of parameters:
            - xi,
            - eta_b, baseline of carbon intensity
            - beta, link between transmissivity and mobility
            - tau_u, UK recovery rate
            - tau_a, adaptation rate
            - tau_l, the rate of lockdown
            - mu_star, critical mu 
            - mu_r, mu recovery
            - alpha, freefall scale factor
            - R_crit, the trigger rate of reproduction for lockdown.
            - R_easy, the trigger rate of reproduction for easing lockdown.
            - sigma, noise strength
        
        - Y0, the starting GDP (units: £billiion)
        
        - Are we in lockdown case (True or False, default is set to False)
        
    Outputs:
        
        - Y_new, the updated GDP (units: £billiion)
        - eta_new, the updated carbon intensity
        - E_new, the updated emissions
        - R_new, the updated reproduction value
        - m_new, the updated mobility
        - mu_new, the updated decarbonisation 
    
    """

    
    # Parameters
    xi = pars[0]
    eta_b = pars[1]
    beta = pars[2]
    tau_u = pars[3]
    tau_a = pars[4]
    tau_l = pars[5]
    mu_star = pars[6]
    mu_r = pars[7]
    alpha = pars[8]
    R_crit = pars[9]
    R_easy = pars[10]
    sigma = pars[11]
    
    # # Initialising arrays
    # Y = np.zeors(n+1)
    # eta = np.zeors(n+1)
    # E = np.zeors(n+1)
    # R = np.zeors(n+1)
    # m = np.zeors(n+1)
    # mu = np.zeors(n+1)
    
    # Initial conditions
    Y = u[0]
    eta = u[1]
    E = u[2]
    R = u[3]
    m = u[4]
    mu = u[5]
    
    # Solving ODE system over a single time step
    Y_new = Y + dt*mu*Y + np.sqrt(dt)*sigma*Y*np.random.normal(0,1)
    eta_new = eta - dt*xi*(eta - eta_b)
    E_new = E + dt*Y*(mu*eta - xi*(eta - eta_b))
    if Lockdown:
        
        R_new = R - dt*R*(beta*tau_l*m + tau_a)
        m_new = m - dt*tau_l*m
        mu_new =  mu - dt*tau_l*mu_r*alpha
        
    else:
        
        R_new = R + dt*R*(beta*tau_u*(1 - m) - tau_a)
        m_new = m + dt*tau_u*(1 - m)
        
        if Y < Y0*np.exp(mu_star*(t-2020-31*1.0/365.25)):
            
            mu_new = min(mu_r,mu + dt*tau_u*(1/m-1)*mu_r*alpha)
            
        else:
            
            mu_new = mu_star
    
    return(Y_new, eta_new, E_new, R_new, m_new, mu_new)
        

def run_simulations():
    
    """
    Run ODE over time period of outbreak to 2050 and return key parameters and variables.
    """
    
    # Time parameters
    dt = 1.0/365.25
    tstart = 2020 + 31*1/365.25
    tend = 2050
    t = np.arange(tstart, tend+dt, step = dt)
    n = len(t)
    
    # Parameters
    xi = 0.089
    eta_b = 0.11
    beta = 2#
    tau_u = 1.57 
    tau_a = 4#
    tau_l = 20.28
    mu_star = 0.0165
    mu_r = 0.05
    alpha = 10
    R_crit = 1.5#
    R_easy = 0.8#
    sigma = 0
    
    # Ket parameters
    tau_us = [0.5, 0.5, 1.75, 3, 3]
    tau_as = [0.5, 3, 1.75, 0.5, 3]
    num_sim = len(tau_us)
    # Initialising arrays
    Y = np.zeros((n+1,num_sim))
    eta = np.zeros((n+1,num_sim))
    E = np.zeros((n+1,num_sim))
    R = np.zeros((n+1,num_sim))
    m = np.zeros((n+1,num_sim))
    mu = np.zeros((n+1,num_sim))
    
    # Initial conditions
    Y[0,:] = 2090
    eta[0,:] = 0.2
    E[0,:] = 415
    R[0,:] = 4
    m[0,:] = 1
    mu[0,:] = mu_star
    
    Lockdown = np.zeros(num_sim,dtype=bool)
    
    for j in range(len(tau_us)):
        
        for i in range(n):
            
            if R[i,j] > R_crit and Lockdown[j] == False:
                Lockdown[j] = True
            elif R[i,j] < R_easy  and Lockdown[j] == True:
                Lockdown[j] = False
                
            pars = [xi, eta_b, beta, tau_us[j], tau_as[j],tau_l, mu_star, mu_r, alpha, R_crit, R_easy, sigma]
                
            u = [Y[i,j], eta[i,j], E[i,j], R[i,j], m[i,j], mu[i,j]]
            
            Y[i+1,j], eta[i+1,j], E[i+1,j], R[i+1,j], m[i+1,j], mu[i+1,j] = ODE_model(t[i], dt, u, pars, Y[0,j], Lockdown[j])
            
    return tau_us, tau_as, Y, eta, E, R, m, mu