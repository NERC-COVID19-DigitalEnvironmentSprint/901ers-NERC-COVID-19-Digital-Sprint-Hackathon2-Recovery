# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 10:19:03 2020

@author: aargl
"""


import numpy as np


def model_lockdown_m(t,m_l,tau_c):

    """
    Function computes the rate of lockdown on the mobility.
    
    Inputs:
        
        - time
    
        - m_l the lower bound of the mobility
        
        - tau_c, the rate of lockdown of mobility.
        
    Outputs:
        
        - m the current mobility
    
    """
        
    return (1-m_l)*np.exp(-tau_c*t)
    

def model_easing_m(t,m_l,tau_l):
    
    """
    Function computes the rate of easing and change in mobility over time.
    
    Inputs:
        
        - time
        
        - m_l the lower bound of the mobility
        
        - tau_l, the rate of easing of lockdown of the system.
    
    Outputs:
        
        - current mobility.        
    
    """
    
    
    return 1 - (1-m_l)*np.exp(-tau_l*t)
    



def model_R(t,R0,beta,m,tau_c):
    
    """
    Function estimates how the R value varys with mobility.
    
    Inputs:
        
        - t, time since R0. (yrs)
        
        - R0, the basic transmission rate.
        
        - beta, the rate of sucepctiblity from the SIR model.
        
        - m, the transmission rate. 
        
        - tau_a, the rate of adaptation.
        
    Outputs:
        
        - an estimate of R. the transmission rate.
    
    """
    
    return R0*np.exp(-(1-m)*beta-tau_a*t)

    


def model_dgdp(Y,dt,mu,sigma):
    
    """
    Estimate the change in gross domestic product (GDP) using a Brownian motion.
    
    Inputs:
        
        - Y, the current GDP (units: £billiion)
        
        - dt, the current time-step (units:yr)
        
        - mu, the annual GDP growth-rate (units:1/yr)
    
        - sigma, the annual volatility/standard deviation (units:1/yr)
        
    Outputs:
        
        - Y + dY, the updated GDP over the given time-step (£billiion)
    
    """
    dW = np.random.randn()*np.sqrt(dt)
    
    dY = mu*Y*dt + sigma*Y*dW

    return Y + dY


def model_expected_gdp(t,Y0,mu):
    
    """
    Estimate the expected mean gross domestic product (GDP), assuming brownian 
    motion volatility.
    
    Inputs:
    
        - t, the current time (units:yr)
        
        - Y0, the initial GDP at t0 (units: £billiion)
    
        - mu, the annual GDP growth-rate (units:%/yr)
    
    Outputs:
        
        - Y0*np.exp(mu*t), the expected GDP over time. (units: £billiion)
    
    """
    
    return Y0*np.exp(mu*t)
        
def model_expected_eta(t,eta_0,eta_b,xi):
    
    """
    Estimate the change in carbon intensity rate (eta) assuming a exponential 
    decay towards an assumed baseline carbon intensit
    
    Inputs:
    
        - t, the current time (units:yr)
        
        - eta_0, the inital carbon intensity (units: MtCO2e/£billion)
        
        - eta_b, the baseline carbon intensity (units: MtCO2e/£billion)
        
        - mu, the decarbonisation rate of GDP (units:%/yr)
        
    
    Outputs:
        
        - (eta_0-eta_b)*np.exp(-xi*t)+eta_b, the current times carbon intensity (units: MtCO2e/£billion)
    
    """
    
    return (eta_0-eta_b)*np.exp(-xi*t)+eta_b


def model_emissions(Y,eta):
    
    """
    Estimate the emissions by assuming that they are equivalent to the product 
    of the GDP and GDP carbon intensity.
    
    Inputs:
        
        - Y, the GDP. (units: £billiion)
        
        - eta, the carbon intensity. (units: MtCO2e/£billion)
        
    Outputs:
        
        - eta * Y, the corresponding total emissions. (units: MtCO2e) 
        
    
    """
    
    return eta * Y
