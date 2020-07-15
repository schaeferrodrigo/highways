# -*- coding: utf-8 -*-
#==============================================================================
import numpy as np
import parameeters as par
import alpha as al
import tau_star as tau
import sys
import newton_method_vv as new

def A_1(I):
    if I != 0:
        a =  2 * np.pi * I *par.a_1/ np.sinh( np.pi * I/2)
    else:
        a = 4 * par.a_1
    return a

def der_A_1(I):
    if I != 0 :
        der_a = par.a_1 * ( (2*np.pi)/(np.sinh(np.pi * I /2)) - (I * np.pi**2)/(np.sinh(np.pi * I/2) * np.tanh(np.pi * I/2)))
    else:
        der_a = 0
    return der_a

def A_2(I):
    if I != 1:
        a =  2 * np.pi * (I - 1) * par.a_2 /( np.sinh( (I -1 ) * np.pi/2 ) )
    else:
        a = 4 * par.a_2
    return a

def der_A_2( I ):
    if I != 1:
        der_a = par.a_2 * ( (2*np.pi)/(np.sinh(np.pi * (I - 1) /2)) - ((I- 1) * np.pi**2)/(np.sinh(np.pi * (I - 1)/2) * np.tanh(np.pi*(I - 1)/2)))
    else:
        der_a = 0
    return der_a

#------------------------------------------
def xi_hor(I , psi):
    """horizontal crest"""
    xi = (-np.arcsin( al.alpha( I ) * par.mu * np.sin( psi ) ))%( 2 * np.pi )
    return xi

def part_der_xi_psi_hor( I, psi ):
    """derivative of the horizontal crest in relation to variable psi"""
    part = -(al.alpha( I ) * (par.mu) * np.cos( psi ) )/ np.sqrt( 1 - ( al.alpha( I ) * (par.mu) * np.sin( psi ) )**2 )
    return part

#-------------------------------------------------------------------------------
def red_poi_fun( I , theta , tau):
    '''funçao reduzida de Ponicaré'''
    red = A_1(I) * np.cos( theta - I * tau) + A_2(I) * np.cos(theta - (I - 1) * tau)
    return red

def dot_I( I , theta, tau ):
    ''' dI/dt'''
    dot_I = A_1( I ) * np.sin( theta - I * tau ) / ( I - 1 )
    return dot_I

def dot_theta( I , theta , tau ):
    '''dtheta/dtheta'''
    dot_theta = der_A_1( I ) * np.cos( theta - I * tau ) + der_A_2( I ) * np.cos( theta - (I - 1)*tau) + tau *( A_1( I ) * np.sin( theta - I * tau) + A_2( I ) * np.sin( theta - (I - 1) * tau))
    return dot_theta

#-------------------------------------------------------------------------------
def tangency( I ,psi):
    '''condição de tangência'''
    tan = part_der_xi_psi_hor( I , psi) - ((I-1)/I)
    return tan

def root_tangency( I , psi_0 , psi_1 , step) :
    '''retorna psi onde ocorre tangência'''
    while tangency( I , psi_0) * tangency( I , psi_1) > 0:
        psi_1 = psi_1 + step
    psi_sec = psi_0 - tangency(I, psi_0)*( psi_0 - psi_1)/(tangency( I , psi_0 ) - tangency( I , psi_1) )
    while np.abs(tangency( I , psi_sec) )> 1e-6:
        if tangency( I , psi_0 ) * tangency( I , psi_sec ) < 0:
            psi_1 = psi_sec
        else:
            psi_0 = psi_sec
        psi_sec = psi_0 - tangency(I, psi_0)*( psi_0 - psi_1)/(tangency( I , psi_0 ) - tangency( I, psi_1) )
    root = psi_sec
    return root

def psi( I ):
    '''retorna os dois valores de psi onde ocorre tangência para um I fixado'''
    step = 0.01
    #primeiro valor de psi
    psi_0 = 0.
    psi_1 =  step
    psi_A = root_tangency( I , psi_0 , psi_1 , step)%(2 * np.pi)
    #segundo valor de psi :
    psi_3 = psi_A + step
    psi_4 = psi_3 + step
    psi_B = root_tangency( I , psi_3 , psi_4 , step )%(2 * np.pi)
    return [ psi_A , psi_B]

def theta_tau_candidate( I , psi , sign  ):
    i = 0
    step = 0.01
    psi_tan = psi
    tau_0 = 0
    if sign == 'pos':
        tau_1 = step
    else:
        tau_1 = -step
    signal = np.sign( tau_1 )
    while tau.eq_crest( I , psi_tan  + I * tau_0 , tau_0  ) * tau.eq_crest( I , psi_tan  + I * tau_1 , tau_1  ) > 0 or ((psi_tan + tau_0 )%(2*np.pi) > np.pi/2 and (psi_tan + tau_0 )%(2*np.pi) < 3*np.pi/2 ):
        tau_0 , tau_1 = tau_1 , tau_1 + signal * step
    # método da bisseçao
    while  np.abs( tau_1 - tau_0 ) > 1e-6:
        tau_c = ( tau_1 + tau_0 )/2
        if tau.eq_crest( I , psi_tan + I * tau_1 , tau_1) * tau.eq_crest( I , psi_tan + I * tau_c , tau_c) < 0:
            tau_0 = tau_c
        else:
            tau_1 = tau_c
    # método da secante
    tau_sec = tau_1 - tau.eq_crest(I, psi_tan + I* tau_1, tau_1)*( tau_1 - tau_0)/(tau.eq_crest( I , psi_tan + I * tau_1, tau_1 ) - tau.eq_crest( I , psi_tan + I * tau_0 , tau_0) )
    i = 0
    while np.abs( tau.eq_crest(I, psi_tan + I * tau_sec , tau_sec) ) > 1e-12:
        i = i + 1
        if i > 1000:
            sys.exit('nao converge')
        if tau.eq_crest( I , psi_tan + I * tau_1 , tau_1 ) * tau.eq_crest( I , psi_tan + I * tau_sec , tau_sec) < 0 :
            tau_0  = tau_sec
        else:
            tau_1 = tau_sec
        tau_sec =  tau_1 - tau.eq_crest(I, psi_tan + I* tau_1 , tau_1)*( tau_1 - tau_0)/(tau.eq_crest( I , psi_tan + I * tau_1 , tau_1 ) - tau.eq_crest( I , psi_tan + I * tau_0 , tau_0) )
    if np.linalg.norm( new.F(I , psi_tan + I*tau_sec, tau_sec ))>1e-4:
        sys.exit('nao eh tangente')
    return  [tau_sec, (psi_tan + I*tau_sec)%(2*np.pi)]

def theta( I , psi ):
    tau_theta_pos = theta_tau_candidate( I , psi , 'pos')
    tau_theta_neg = theta_tau_candidate( I , psi , 'neg' )
    tau_pos = tau_theta_pos[0]
    tau_neg = tau_theta_neg[0]
    if np.minimum( np.abs( tau_pos) , np.abs( tau_neg )) == np.abs(tau_pos):
        tau_final, theta_final = tau_theta_pos[0] , tau_theta_pos[1]
    else:
        tau_final, theta_final = tau_theta_neg[0] , tau_theta_neg[1]
    if np.linalg.norm( new.F(I , theta_final, tau_final))>1e-4:
        print np.linalg.norm( new.F(I , psi  + I*tau_neg, tau_neg))>1e-4 , np.linalg.norm( new.F(I , psi  + I*tau_neg, tau_pos))>1e-4
        sys.exit('nao eh tangente')
    return [theta_final , tau_final]
