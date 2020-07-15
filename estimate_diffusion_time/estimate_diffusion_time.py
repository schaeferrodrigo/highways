# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
from alpha_function import *
from scipy import special,optimize


def diffusion_time(T_s , angle , A , eps):
    "estimate of the diffusion_time time. ref=regular and chaotic dynnmics"
    cstnt = 16 * np.cos(angle) *( 1 + 1.465/np.sqrt(1-(A/np.tan(angle))**2)) #constante C
    dif_time = (T_s/eps) * 2 * np.log(cstnt/eps)
    return dif_time
