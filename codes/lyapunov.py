#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 01:37:47 2021

@author: nehabinish
"""


import numpy as np
import matplotlib.pyplot as plt
import numba


#DEFINING THE RUNGE KUTTA 4 INTEGRATOR

def RK4(f, x0, t):
    '''
    

    Parameters
    ----------
    f : function to be integrated 
    
    x0 : initial conditions 
    
    t0 : initial time 
    
    t1 : final time 
    
    dt : time step 

    Returns
    -------
    X : An array of derivatives 

    '''
    
    dt = t[1] - t[0]
    N = len(t)
    X = np.empty((len(t), len(x0)))
    X[0] = x0
    
    for i in range(1, N):
        
        k1 = f(X[i-1], t[i-1])
        k2 = f(X[i-1] + dt/2*k1, t[i-1] + dt/2)
        k3 = f(X[i-1] + dt/2*k2, t[i-1] + dt/2)
        k4 = f(X[i-1] + dt*k3, t[i-1] + dt)
        
        X[i] = X[i-1] + dt/6*(k1 + 2*k2 + 2*k3 + k4)
        
    return X
 


# DEFINING THE CHAOTIC SEIR MODEL AS A FUNCTION 

def func(y, t):
    '''
    

    Parameters
    ----------
    t : time vector 
    
    y : an array with all the compartments 
    
    beta : parameter beta
    
    sigma : parameter sigma 
    
    gamma : parameter gamma 

    Returns
    -------
    TYPE
        array ( ds/dt, de/dt, di/dt, dr/dt )

    '''

    
    S, E, I, R = y

    
    beta  = beta0 * ( 1 + (beta1 * np.cos(2*np.pi*t)))
    
    dsdt  = b - (b*S)  - (beta*S*I)
    dedt  = (beta*S*I) - ( (alpha+b)*E )
    didt  = (alpha*E)- ((gamma +b)*I)
    drdt  = (gamma * I) - (b*R)

    
    return np.array([ dsdt, dedt, didt, drdt ])

#Jacobian matrix
def JM(v, t):
    
    S, E, I, R = [a for a in v]
    
    beta  = beta0 * ( 1 + (beta1 * np.cos(2*np.pi*t)))
        
    Df = [ [ -(b+ (beta * I)), 0               , - beta * S     , 0  ], 
           [ beta * I        , - (alpha + beta),   beta * S     , 0  ],  
           [ 0               , alpha           , -(gamma + beta), 0  ],  
           [ 0               , 0               ,   gamma        , -b ] ] 
    
    return np.array(Df)




beta1 = 0.9
b     = 0.02
alpha = 35.84
gamma = 100
beta0 = 1800
    
s0 = 0.6
e0 = 0.3
i0 = 0.05
r0 = 0.05
    
        
x0 = [ s0, e0, i0, r0 ]
t0 = [0,100]
        
dt = 0.001    
t = np.arange(t0[0], t0[1], dt)   #time span 
N = len(t)
    
v_n = RK4(func, x0, t)


U = np.eye(4)  

lyap = [] # list for lyapunov

# QM decomposition
for k in range(0, N):
    
    v0 = v_n[k] # updating v0 after every iteration

    U_n = np.matmul( np.eye(4) + JM(v0, t[k]) * dt, U)

    # Gram-Schmidt Orthogonalisation
    
    Q, R = np.linalg.qr(U_n)
    
    lyap.append(np.log(abs(R.diagonal())))

    U = Q #new axes after iteration
    
    print(lyap[k])

print([sum([lyap[l][j] for l in range(N)]) / (dt * N) for j in range(4)])









