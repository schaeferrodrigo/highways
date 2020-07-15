# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
from alpha_function import *


def eq_crest(I , theta , tau, angle):
    "crest's equation"
    crista = alpha_function(I)*(1/np.tan(angle)) * np.sin( theta - I * tau) + np.sin( - tau)
    return crista

def bissec_method( I , theta , tau_1 , tau_2 , angle ,tol = 1e-2):
    "bissec method"
    while np.abs( tau_1 - tau_2 ) > tol:
        tau_c = ( tau_1 + tau_2 )/2
        if eq_crest( I , theta , tau_1,angle) * eq_crest( I , theta , tau_c,angle) < 0:
            tau_2 = tau_c
        else:
            tau_1 = tau_c
    return [tau_1 , tau_2]

def secant_method(I , theta, tau_1 , tau_2 , angle , tol = 1e-10):
    tau_sec = tau_1 - eq_crest(I, theta, tau_1 , angle)*( tau_1 - tau_2)/(eq_crest( I , theta, tau_1 ,angle ) - eq_crest( I , theta , tau_2, angle) )
    while np.abs( eq_crest(I, theta , tau_sec , angle ) ) > tol:
        if eq_crest( I , theta , tau_1 ,angle ) * eq_crest( I , theta , tau_sec, angle) < 0 :
            tau_2  = tau_sec
        else:
            tau_1 = tau_sec
        tau_sec =  tau_1 - eq_crest(I, theta, tau_1 , angle )*( tau_1 - tau_2)/(eq_crest( I , theta, tau_1 , angle ) - eq_crest( I , theta , tau_2 , angle) )
    return tau_sec

def tau(I , theta ,sign , angle , tau_initial = 0 , step = 0.01):
    tau_2 = tau_initial
    if sign == 'neg':
        tau_1 = -step
    else:
        tau_1 = step
    sign = np.sign(tau_1)
    while eq_crest( I , theta , tau_1 , angle) * eq_crest( I , theta , tau_2 ,angle) > 0 :
        tau_2 = tau_1
        tau_1 = tau_1 + sign * step
    candidate = bissec_method( I , theta , tau_1 , tau_2 ,angle)
    tau_1 = candidate[0]
    tau_2 = candidate[1]
    tau_star = secant_method(I , theta, tau_1 , tau_2,angle) #método secanda
    return tau_star

def assign_tau( I, theta ,angle ):
    value_of_tau_pos = tau( I , theta , 'pos' , angle )
    value_of_tau_neg = tau( I , theta , 'neg' , angle )
    #print value_of_tau_neg , value_of_tau_pos
    #print 'value', I
    #print 'candidatos a tau' , value_of_tau_neg , value_of_tau_pos
    if np.minimum( np.abs(value_of_tau_neg) , np.abs(value_of_tau_pos) ) == np.abs(value_of_tau_neg):
        val_of_tau = value_of_tau_neg
    else:
        val_of_tau = value_of_tau_pos
    return val_of_tau
