# -*- coding: utf-8 -*-
#==============================================================================
# Scattering map on I and Theta
#==============================================================================
import numpy as np

def alpha( I ):
    '''definiçao da funçao alpha'''
    if I != 0 :
        alpha = (np.sinh( (I-1) * np.pi/2 )/np.sinh( np.pi * I /2)) * (I/(I-1))**2
    else:
        alpha = 0
    return alpha

def beta(I):
    '''definiçao da funçao beta'''
    beta = I * alpha(I)/(I-1)
    return beta
