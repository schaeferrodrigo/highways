# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np

def A_1(I , angle):
    a_1 = 2*np.pi*np.cos(angle)*I /np.sinh(np.pi * I/2)
    return a_1

def A_2(angle):
    a_2 = 2* np.pi*np.sin(angle)/np.sinh(np.pi/2)
    return a_2

def func(I , angle):
    fun = ((I**2)*A_2(angle) - np.sqrt(A_2(angle)**2 + (I**2 - 1) * (I**2)*A_1(I,angle)**2))/((I**2 - 1)*A_2(angle))
    #print 'fun = ', 1- fun
    return fun

def theta_higways( I , angle):
    "theta as a function of I"
    print 'teste= ', (A_2(angle)/A_1(I , angle)) * (1-func(I,angle))
    theta = - np.arccos((A_2(angle)*(1 - func(I,angle)))/A_1(I , angle)) + I * np.arccos(func(I,angle))
    return theta
