#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 19:30:16 2021

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
def seir_rand(y, t):
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


beta_val = np.linspace(0.1,0.9,5)

for beta1 in beta_val:

    # parameters

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


    t = np.arange(t0[0], t0[1], 0.001)   #time span
    #t = np.linspace(0,10,500)


    parameters = RK4(seir_rand, x0, t)



    fig = plt.figure(figsize = (10,5));

    plt.plot(parameters[20000:,0], parameters[20000:,1], '-')
    plt.xlabel('S(t)')
    plt.ylabel('E(t)')
    plt.title(' beta1 = {}'.format(beta1))
    plt.grid()


    fig = plt.figure(figsize = (10,5));

    plt.plot(parameters[20000:,0], parameters[20000:,2], '-')
    plt.xlabel('S(t)')
    plt.ylabel('I(t)')
    plt.title(' beta1 = {}'.format(beta1))
    plt.grid()

    fig = plt.figure(figsize = (10,5));

    plt.plot(parameters[20000:,1], parameters[20000:,2], '-')
    plt.xlabel('E(t)')
    plt.ylabel('I(t)')
    plt.title(' beta1 = {}'.format(beta1))
    plt.grid()



    fig = plt.figure(figsize = (10,5));
    ax = plt.axes(projection='3d')

    # Data for a three-dimensional line
    zline = parameters[20000:,0]
    xline = parameters[20000:,1]
    yline = parameters[20000:,2]
    ax.plot3D(xline, yline, zline, 'red')
    ax.set_xlabel('S')
    ax.set_ylabel('E')
    ax.set_zlabel('I');
    ax.set_title(' beta1 = {}'.format(beta1))


    #fig = plt.figure(figsize = (10,5));
    #ax = plt.axes(projection='3d')

    # Data for a three-dimensional line
    #zline = parameters[:,0]
    #xline = parameters[:,2]
    #yline = parameters[:,3]
    #ax.plot3D(xline, yline, zline, 'red')
    #ax.set_xlabel('S')
    #ax.set_ylabel('I')
    #ax.set_zlabel('R');
    #ax.set_title(' beta1 = {}'.format(beta1))

    #fig = plt.figure(figsize = (10,5));
    #ax = plt.axes(projection='3d')

    # Data for a three-dimensional line
    #zline = parameters[:,0]
    #xline = parameters[:,1]
    #yline = parameters[:,3]
    #ax.plot3D(xline, yline, zline, 'red')
    #ax.set_xlabel('S')
    #ax.set_ylabel('E')
    #ax.set_zlabel('R')

    #ax.set_xlim(0,0.3)
    #ax.set_ylim(0,0.3)

    #ax.set_title(' beta1 = {}'.format(beta1))

    #lyapunov1.append(np.mean(result1))
    #lyapunov2.append(np.mean(result2))
    #lyapunov3.append(np.mean(result3))
    #lyapunov4.append(np.mean(result4))

    fig = plt.figure(figsize = (10,5)); ax = fig.gca()
    ax.plot(t[20000:], parameters[20000:,0], '-', label = 'S')
    ax.plot(t[20000:], parameters[20000:,1], '-', label = 'E')
    ax.plot(t[20000:], parameters[20000:,2], '-', label = 'I')
    ax.plot(t[20000:], parameters[20000:,3], '-', label = 'R')
    plt.grid()
    plt.legend()
    plt.xlabel('Time - Number of days ')
    plt.ylabel('Number of people')
    plt.title( 'beta1 = {}'.format(beta1))


    fig = plt.figure(figsize = (10,5)); ax = fig.gca()
    ax.plot(t, parameters[:,0], '-', label = 'S')
    ax.plot(t, parameters[:,2], '-', label = 'I')
    plt.grid()
    plt.legend()
    plt.xlabel('Time - Number of days ')
    plt.ylabel('Number of people')
    plt.title( 'beta1 = {}'.format(beta1))
