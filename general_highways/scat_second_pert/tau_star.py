# -*- coding: utf-8 -*-
#==============================================================================
import numpy as np
import functions as fun

def eq_crest(I , theta , tau):
    crista = I * fun.A_1( I ) * np.sin( theta - I * tau) + ( I - 1 ) * fun.A_2( I ) * np.sin( theta - (I - 1) * tau)
    return crista

def bissec_method( I , theta , tau_1 , tau_2  , tol = 1e-2):
    while np.abs( tau_1 - tau_2 ) > tol:
        tau_c = ( tau_1 + tau_2 )/2
        if eq_crest( I , theta , tau_1) * eq_crest( I , theta , tau_c) < 0:
            tau_2 = tau_c
        else:
            tau_1 = tau_c
    return [tau_1 , tau_2]

def secant_method(I , theta, tau_1 , tau_2 , tol = 1e-10):
    tau_sec = tau_1 - eq_crest(I, theta, tau_1)*( tau_1 - tau_2)/(eq_crest( I , theta, tau_1 ) - eq_crest( I , theta , tau_2) )
    while np.abs( eq_crest(I, theta , tau_sec) ) > tol:
        if eq_crest( I , theta , tau_1 ) * eq_crest( I , theta , tau_sec) < 0 :
            tau_2  = tau_sec
        else:
            tau_1 = tau_sec
        tau_sec =  tau_1 - eq_crest(I, theta, tau_1)*( tau_1 - tau_2)/(eq_crest( I , theta, tau_1 ) - eq_crest( I , theta , tau_2) )
    return tau_sec

def tau(I , theta ,sign , tau_initial = 0 , step = 0.05):
    tau_2 = tau_initial
    if sign == 'neg':
        tau_1 = -step
    else:
        tau_1 = step
    sign = np.sign(tau_1)
    while eq_crest( I , theta , tau_1 ) * eq_crest( I , theta , tau_2 ) > 0 :
        tau_2 = tau_1
        tau_1 = tau_1 + sign * step
    candidate = bissec_method( I , theta , tau_1 , tau_2 )
    tau_1 = candidate[0]
    tau_2 = candidate[1]
    tau_star = secant_method(I , theta, tau_1 , tau_2) #método secanda
    return tau_star
