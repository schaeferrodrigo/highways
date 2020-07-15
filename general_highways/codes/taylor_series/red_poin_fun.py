# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
import matplotlib.pyplot as plt
from parametros import *
from functions import *
from find_tau import *
import bigfloat as bf
#==============================================================================
def red_poi_fun( I , theta , tau ,r=0):
    '''funçao reduzida de Poincaré'''
    red = A(I,a_1) * bf.cos( theta - I * tau) + A(I*r - 1,a_2) * bf.cos( r*theta -(I*r-1)* tau)
    return red

#recalcular
def dot_I( I , theta, tau ,r):
    ''' dI/dt'''
    dot_I = - A( I , a_1 ) * bf.sin( theta - I * tau )/(I*r-1)
    return dot_I

def part_red_poin_r(I,theta,tau):
    r_0 = I*der_A(-1,a_2)*bf.cos(tau) - A(-1,a_2)*(theta - I*tau)*bf.sin(tau)
    return r_0

def red_poin_approx(I,theta, tau):
    first_term = A(-1,a_2)
    second_term = r * part_red_poin_r(I,theta,tau)
    sec_part_red_poin_r =( (I**2)* der_2_A(-1,a_2)*bf.cos(tau)
                         -der_A(-1,a_2)*I*(theta-I*tau)*bf.sin(tau)
                         -der_A(-1,a_2)*tau_dr(I,theta,tau)*bf.sin(tau)
                         -I*der_A(-1,a_2)*(theta-I*tau)*bf.sin(tau)
                         +A(-1,a_2)*I*tau_dr(I,theta,tau)*bf.sin(tau)
                         -A(-1,a_2)*((theta - I*tau)**2)*bf.cos(tau)
                         -A(-1,a_2)*(theta-I*tau)*tau_dr(I,theta,tau)*bf.cos(tau))
    third_term = (r**2)* sec_part_red_poin_r
    return first_term + second_term + third_term
