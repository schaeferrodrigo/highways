# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
import matplotlib.pyplot as plt
from parametros import *
from functions import *
#===============================================================================
 # método numerico encontrar tau

def eq_crest(I , theta , tau):
    crista = I * A( I , a_1 ) * np.sin( theta - I * tau) + (I*r-1) *A( I*r -1 , a_2 ) * np.sin(r*theta- (I*r-1) *tau)
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
    tau_sec = tau_1 - eq_crest(I, theta, tau_1)*( tau_1 - tau_2)/(eq_crest( I , theta, tau_1) - eq_crest( I , theta , tau_2) )
    while np.abs( eq_crest(I, theta , tau_sec) ) > tol:
        if eq_crest( I , theta , tau_1 ) * eq_crest( I , theta , tau_sec) < 0 :
            tau_2  = tau_sec
        else:
            tau_1 = tau_sec
        tau_sec =  tau_1 - eq_crest(I, theta, tau_1)*( tau_1 - tau_2)/(eq_crest( I , theta, tau_1) - eq_crest( I , theta , tau_2) )
    return tau_sec

def tau(I , theta ,sign , tau_initial = 0 , step = 0.05):
    tau_2 = tau_initial
    if sign == 'neg':
        tau_1 = -step
    else:
        tau_1 = step
    sign = np.sign(tau_1)
    altura = 0
    while (eq_crest( I , theta , tau_1 ) * eq_crest( I , theta , tau_2 ) )> 0 :
            tau_2 = tau_1
            tau_1 = tau_1 + sign * step
    candidate = bissec_method( I , theta , tau_1 , tau_2 )
    cand_0 = candidate[0]
    cand_1 = candidate[1]
    tau_star_cand = secant_method(I , theta, cand_0 , cand_1)
    altura =  -tau_star_cand
    tau_star = tau_star_cand
    return tau_star

def assign_tau( I, theta ):
    tau_pos = tau( I , theta , 'pos')
    tau_neg = tau( I , theta , 'neg' )
    if np.minimum( np.abs(tau_neg) , np.abs(tau_pos) ) == np.abs(tau_neg):
        val_of_tau = tau_neg
    else:
        val_of_tau = tau_pos
#     #if np.abs( val_of_tau ) > np.pi/2:
#     #     sys.exit(' estou intersectando a crista errada')
    return val_of_tau
#==============================================================================
