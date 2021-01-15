# -*- coding: utf-8 -*-
#===============================================================================
from parametros import *
import numpy as np
from find_tau import *
#==============================================================================

def function_to_be_bounded(theta):
    tau_star = assign_tau(I ,theta)
    return np.abs((r*I-1)*np.sqrt(1-(mu*alpha()*np.sin(theta - I*tau_star))**2)+I*mu*alpha()*np.cos(theta-I*tau_star))
