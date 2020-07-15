# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
import matplotlib.pyplot as plt
from parametros import *
#===============================================================================
# Funcoes

def A(I,a):
    """
    """
    A = 2*np.pi*I*a/np.sinh(np.pi*I/2)
    return A

def der_A(I,a):
    """
    """
    der =( 2 * np.pi* a * ( np.sinh( np.pi *I / 2) -
     I * (np.pi / 2)* np.cosh( np.pi *I / 2))/(np.sinh(np.pi * I /2)**2))
    return der

def der_2_A(I,a):
    der_2_num = (-(I*np.pi**2/4)*np.sinh(np.pi*I/4)**2
                - np.sinh(np.pi*I/2)*np.cosh(np.pi*I/2)*np.pi
                +(I*np.pi**2/2)*np.cosh(np.pi*I/2)**2)
    der_2_den = np.sinh(np.pi*I/2)**3
    return 2*np.pi*a*der_2_num/der_2_den            


def tau_dr( I , theta , tau):
    num = ((-A(I,a_1) + der_A(-1,a_2))*I*np.sin(tau)
          + (theta-I*tau)*A(-1, a_2)*np.cos(tau))
    den = -A(-1,a_2)*(I**2 + (1-I**2)*np.cos(tau))
    return num/den
