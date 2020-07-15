# -*- coding: utf-8 -*-
#=============================================================================
import numpy as np
import parameeters as par
import alpha as al

def sing_value( I , mu ):
    sing = np.abs(mu * al.alpha(I) ) - 1
    return sing

def bissec_method( I_0 , I_1 , mu , tol = 1e-2):
    while np.abs( I_0 - I_1) > tol:
        I_c = (I_0  + I_1 ) / 2
        if sing_value( I_0 , mu ) * sing_value(I_c , mu) < 0:
            I_1 = I_c
        else:
            I_0 = I_c
    return [ I_0 , I_1 ]

def secant_method( I_0 , I_1 , mu , tol = 1e-10 ):
    I_sec = I_1 - sing_value( I_1 , mu )*(I_1 - I_0)/ ( sing_value( I_1 , mu) - sing_value(I_0 , mu))
    while sing_value( I_sec , mu) > tol:
        if sing_value( I_1  , mu ) * sing_value( I_sec , mu ) < 0:
            I_0 = I_sec
        else:
            I_1 = I_sec
        I_sec = I_1 - sing_value( I_1 , mu )*(I_1 - I_0)/ ( sing_value( I_1 , mu) - sing_value(I_0 , mu))
    return I_sec

def singularity_in_I( domain_I , mu , step = 0.01 ):
    if np.abs( mu ) <= np.exp( - np.pi/2) or np.abs( mu ) >= np.exp( np.pi/2):
        cardinalidade = 2
    else:
        cardinalidade = 3
    singul_values_I  = []
    I_0 , I_1  = domain_I[ 0 ] , domain_I[0] + step
    while len( singul_values_I ) != cardinalidade:
        while sing_value( I_0 , mu ) * sing_value(I_1 , mu) > 0:
            I_0  , I_1 =  I_1 , I_1 + step
        candidates_I = bissec_method( I_0 , I_1 , mu)
        I_A , I_B = candidates_I[ 0 ] , candidates_I[ 1 ]
        I_root = secant_method( I_A , I_B , mu )
        singul_values_I.append(I_root)
        if I_1 < 0:
            I_0 , I_1 = 0 , step
        else:
            I_0  , I_1 = I_1 , I_1 + step
    return singul_values_I
