# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
from decimal import *
import  matplotlib.pyplot as plt

gamma = 0.1# test
I_init = 1e-5

def bissec_method(f,tol = 1e-4):
    I_neg, I_pos = 1e-2 , 2
    while f(I_neg,gamma)*f(I_pos,gamma)>0:
        I_neg , I_pos = I_neg-1e-4 , I_pos + 1e-2
    I_1 , I_2 = I_neg , I_pos
    while np.abs(I_1-I_2 ) > tol:
        I_c = (I_1 + I_2)/2
        if f(I_1,gamma)*f(I_c,gamma)<0:
            I_2 = I_c
        else:
            I_1 = I_c
    return [I_1, I_2]

def secant_method(f, I_1 , I_2, gamma , tol = 1e-10):
    I_sec = I_1 - f(I_1 ,gamma)*(I_1 - I_2)/(f(I_1 ,gamma) - f(I_2,gamma))
    while np.abs(f(I_sec,gamma))>tol:
        if f(I_1,gamma) * f(I_sec,gamma) < 0:
            I_2 = I_sec
        else:
            I_1 = I_sec
        I_sec = I_1 - f(I_1, gamma)*(I_1 - I_2)/(f(I_1,gamma) - f(I_2,gamma))
    return I_sec

def to_find_I_plus(f,gamma):
    I_1 , I_2 = bissec_method(f)
    print("bissec finished")
    I_plus = secant_method(f , I_1 , I_2 , gamma)
    print("secant finished")
    return I_plus


def fun_I_plus_1(I , gamma):
    gamma = np.float128(gamma)
    return ((I**3)*np.sinh(np.pi/2)/np.sinh(np.pi*I/2)) - 1/np.abs(np.tan(gamma*np.pi/2))

def fun_I_plus_2(I , gamma):
    gamma = np.float128(gamma)
    return ((I**2)*np.sinh(np.pi/2)/np.sinh(np.pi*I/2)) - 1/np.abs(np.tan(gamma*np.pi/2))


def domain_I(I_init, gamma,step = 0.01 ):
    mu_abs = np.abs(np.tan(np.float128(gamma*np.pi/2)))
    print("mu = " , mu_abs)
    if mu_abs< 0.625:
        I_plus  = 5
    else:
        if mu_abs <=1:
            I_plus = to_find_I_plus(fun_I_plus_1,gamma)
        else:
            I_plus = to_find_I_plus(fun_I_plus_2 , gamma)
    return  np.arange(float(I_init) , I_plus , step)

def alpha(I):
    "alpha function that comes from crests' equation"
    #I = Decimal(str(I))
    return ((I**2) * np.sinh(np.float128(np.pi/2)))/np.sinh(np.float128(I*np.pi/2))

def function_f( I , gamma):
    "Function f defined in order to simply the expression of \tau^* function"
    mu = np.tan(np.float128(gamma*np.pi/2))
    return (((I**2) - np.sqrt((1 + (I**2 -1)*(mu**2)*(alpha(I)**2))))/(I**2 -1))



def tau(I,gamma):
    return np.abs(np.arccos(np.float128(function_f(I,gamma))))

def der_func_f( I , gamma):
    num =(I*((np.tan(gamma*np.pi/2)**2)*(I**2)*(I**2 - 1)*(np.sinh(np.pi/2)**2)*(-2*(I**2) + np.pi*(I**2 - 1)*I*(1/np.tanh((np.pi*I)/2)) + 4)* (1/np.sinh((np.pi*I)/2)**2) - 4* np.sqrt((np.tan(gamma*np.pi/2)**2)*(I**2 - 1)*(I**4)*(np.sinh(np.pi/2)**2)*(1/np.sinh((np.pi*I)/2)**2) + 1) + 4))
    den = (2*((I**2 - 1)**2)*np.sqrt((np.tan(gamma*np.pi/2)**2)*(I**2 - 1)*(I**4)*(np.sinh(np.pi/2)**2)*(1/np.sinh((np.pi*I)/2)**2) + 1))
    return (num/den)



def abs_der_fun(I,gamma):
    return (np.abs(der_func_f(I,gamma)/(np.sqrt(1-function_f(I,gamma)**2))))



def I_der_tau(I, gamma):
    return I * abs_der_fun(I,gamma)



I_values = domain_I(I_init,gamma)
plt.plot(I_values, function_f(I_values,gamma), color = 'black' , label=r'$f$')
plt.plot(I_values, tau(I_values,gamma), color= 'red', label = r'$\tau$')
plt.plot(I_values , abs_der_fun(I_values,gamma) ,color = 'blue', label = r"$d\tau/dI$")
plt.plot(I_values , I_der_tau(I_values,gamma) , color = 'green' , label = r"$Id\tau/dI$")
plt.legend()
plt.show()
