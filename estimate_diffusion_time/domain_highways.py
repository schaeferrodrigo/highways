# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
from alpha_function import *
from scipy import special,optimize

def def_interval_I(angle):
    "definition of the interval where the highways are defined ref: regular and chaotic dynamics "
    mu = 1/np.tan(angle)
    if mu < 0.625:
        interval_I = np.array([0.5, 1 , 1.5 , 2])
    elif mu >= 0.625 and mu <= 1.:
        def control_1(I):
            return I *alpha_function(I) - np.tan(angle)
        I_max = optimize.brentq(control_1, a= 0.98 , b= 1.9)
        if I_max < 1.:
            interval_I = np.array([0.5 , 0, 0 ,0])
        elif I_max >=1. and I_max <1.5:
            interval_I = np.array([0.5 , 1. , 0 ,0])
        else:
            interval_I = np.array([0.5 , 1., 1.5 ,0])
    else:
        def control_2(I):
            return alpha_function(I) - np.tan(angle)
        I_max =  optimize.brentq( control_2, a= 0.001 , b= 1.1)
        if I_max >=0.5:
            interval_I = np.array([0.5 , 0 , 0, 0])
        else:
            interval_I = np.array([0]*4)
    return interval_I
