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
            - mu_star, critical mu 
            - mu_r, mu recovery
            - mu_f, mu freefall
            - R_crit
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
    tau_l = pars[4]
    tau_a = pars[5]
    mu_star = pars[6]
    mu_r = pars[7]
    mu_f = pars[8]
    sigma = pars[9]
    
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
        
        if Y < Y0*np.exp(mu_star*(t-2020)):            
            mu_new = mu - dt*tau_l*m*mu_r/(m + mu_f)            
        else:            
            mu_new = mu
            
    else:
        
        R_new = R + dt*R*(beta*tau_u*(1 - m) - tau_a)
        m_new = m + dt*tau_u*(1 - m)
        
        if Y < Y0*np.exp(mu_star*(t-2020)):            
            mu_new = mu + dt*tau_u*(1 - m)*mu_r/(m + mu_f)            
        else:            
            mu_new = mu
                
    return(Y_new, eta_new, E_new, R_new, m_new, mu_new)
        


# Time parameters
dt = 1/365.25
tstart = 2020 + 31*dt
tend = 2050
t = np.arange(tstart, tend+dt, step = dt)
n = len(t)

# Parameters
xi = 0.089
eta_b = 0.12 
beta = 0.1#
tau_u = 1.57 
tau_l = 20#
tau_a = 4#
mu_star = 0.0165
mu_r = 1#
mu_f = 0.3#
R_crit = 2#
R_easy = 0.7#
sigma = 0

# Initialising arrays
Y = np.zeros(n+1)
eta = np.zeros(n+1)
E = np.zeros(n+1)
R = np.zeros(n+1)
m = np.zeros(n+1)
mu = np.zeros(n+1)

# Initial conditions
Y[0] = 2029.53
eta[0] = 0.22
E[0] = 438.83
R[0] = 4
m[0] = 1
mu[0] = mu_star

tau_us = [0.5, 0.5, 1.5, 3, 5]
tau_as = [0.5, 3, 1.5, 0.5, 3]

plt.figure(2)
plt.plot(t, Y[0]*np.exp(mu_star*(t-2020)),'k')

for j in range(len(tau_us)):
    for i in range(n):
        
        Lockdown = False
        
        if R[i] > R_crit and not Lockdown:
            Lockdown = True
        elif R[i] < R_easy and Lockdown:
            Lockdown = False
            
        pars = [xi, eta_b, beta, tau_us[j], tau_l, tau_as[j], mu_star, mu_r, mu_f, sigma]
            
        u = [Y[i], eta[i], E[i], R[i], m[i], mu[i]]
        
        Y[i+1], eta[i+1], E[i+1], R[i+1], m[i+1], mu[i+1] = ODE_model(t[i], dt, u, pars, Y[0], Lockdown)
    
    plt.figure(1)
    plt.plot(t, E[:-1], label = '$\tau_u = $'+str(tau_us[j])+', $\tau_a = $'+str(tau_as[j]))
    plt.figure(2)
    plt.plot(t, Y[:-1], label = '$\tau_u = $'+str(tau_us[j])+', $\tau_a = $'+str(tau_as[j]))
    plt.figure(3)
    plt.plot(t, mu[:-1])
    plt.figure(4)
    plt.plot(t, m[:-1])
    plt.figure(5)
    plt.plot(t, R[:-1])
        
   