# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
import matplotlib.pyplot as plt
from parametros import *
from functions import *
from find_tau import *
#==============================================================================
def red_poi_fun( I , theta , tau):
    '''funçao reduzida de Poincaré'''
    red = A(I,a_1) * np.cos( theta - I * tau) + A(I*r - 1,a_2) * np.cos( r*theta -(I*r-1)* tau)
    return red

#recalcular
def dot_I( I , theta, tau ):
    ''' dI/dt'''
    dot_I = - A( I , a_1 ) * np.sin( theta - I * tau )/(I*r-1)
    return dot_I

def part_red_poin_r(I,theta,tau):
    r_0 = I*der_A(-1,a_2)*np.cos(tau) - A(-1,a_2)*(theta - I*tau)*np.sin(tau)
    return r_0

def red_poin_approx(I,theta, tau):
    first_term = A(-1,a_2)
    second_term = r * part_red_poin_r(I,theta,tau)
    sec_part_red_poin_r =( (I**2)* der_2_A(-1,a_2)*np.cos(tau)
                         -der_A(-1,a_2)*I*(theta-I*tau)*np.sin(tau)
                         -der_A(-1,a_2)*tau_dr(I,theta,tau)*np.sin(tau)
                         -I*der_A(-1,a_2)*(theta-I*tau)*np.sin(tau)
                         +A(-1,a_2)*I*tau_dr(I,theta,tau)*np.sin(tau)
                         -A(-1,a_2)*((theta - I*tau)**2)*np.cos(tau)
                         -A(-1,a_2)*(theta-I*tau)*tau_dr(I,theta,tau)*np.cos(tau))
    third_term = (r**2)* sec_part_red_poin_r
    return first_term + second_term + third_term
