# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
import matplotlib.pyplot as plt
from parametros import *
from functions import *
from find_tau import *
#==============================================================================
def red_poi_fun( I , theta , tau):
    '''funçao reduzida de Ponicaré'''
    red = A(I,a_1) * np.cos( theta - I * tau) + A(1.,1.) * np.cos(  - tau)
    return red

def dot_I( I , theta, tau ):
    ''' dI/dt'''
    dot_I = - A( I , a_1 ) * np.sin( theta - I * tau )
    return dot_I
