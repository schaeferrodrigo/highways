# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
from theta_highways import *
from tau_star import *

def integrand(I,angle):
    theta = theta_higways(I,angle)
    #print 'theta= ' ,theta%(2*np.pi)
    tau_star = assign_tau( I, theta ,angle )
    #print 'psi = ', theta - I*tau_star
    int = -np.sinh(np.pi * I /2)/(np.pi*I*np.sin(angle)*np.sin(theta - I*tau_star))
    return int
