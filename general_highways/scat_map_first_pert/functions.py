# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
import matplotlib.pyplot as plt
from parametros import *
#===============================================================================
# Funcoes

def A(I,a):
    A = 2*np.pi*I* a/np.sinh(np.pi*I/2)
    return A

def der_A(I,a):
    der = 2 * np.pi * a * ( np.sinh( np.pi *I / 2) - I * (np.pi / 2)* np.cosh( np.pi *I / 2)) /(np.sinh(np.pi * I /2)**2)
    return der
    
