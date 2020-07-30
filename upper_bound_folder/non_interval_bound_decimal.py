# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
import sys
sys.path.append("/usr/local/lib/python3.8/site-packages")
from interval import interval,inf,imath
from decimal import *

gamma = Decimal('0.5')# test
I_init = Decimal(str(1e-11))

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
    mu_abs = np.abs(np.tan(np.float128(gamma*Decimal(np.pi)/Decimal('2'))))
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
    I = Decimal(str(I))
    return ((I**2) * Decimal(str(np.sinh(np.float128(Decimal(np.pi/2))))))/Decimal(str(np.sinh(np.float128(I*Decimal(np.pi)/2))))

def function_f( I , gamma):
    "Function f defined in order to simply the expression of \tau^* function"
    I = Decimal(str(I))
    if I != 1 :
        mu = Decimal(str(np.tan(np.float128(gamma*Decimal(np.pi)/2))))
        return (((I**2) - Decimal((1 + (I**2 -1)*(mu**2)*(alpha(I)**2))).sqrt())/(I**2 -1))
    else:
        return 1 - (np.tan(gamma*np.pi/2)**2)*(alpha(I)**2)/2




def max_value(f , I_init, gamma):
    vmax_candidate = f(I_init,gamma)
    pmax_candidate = I_init
    I_values = domain_I(I_init,gamma)
    for I in I_values:
        image_I = f(I,gamma)
        print("value of I = ",I," image= " , image_I)
        if vmax_candidate < image_I:
            print("max_candidate = ", vmax_candidate)
            vmax_candidate = image_I
            pmax_candidate = I
        else:
            pass
    return [vmax_candidate, pmax_candidate]



#print(max_value(function_f,I_init,gamma))


def der_func_f( I , gamma):
    if I != 0 and I != 1:
        mu = Decimal(str(np.tan(np.float128(gamma*Decimal(np.pi)/2))))
        I = Decimal(str(I))
        num =(I*((mu**2)*(I**2)*(I**2 - 1)*(Decimal(str(np.sinh(np.pi/2)))**2)*(-2*(I**2) + Decimal(np.pi)*(I**2 - 1)*I*(1/Decimal(str(np.tanh(np.float128((Decimal(np.pi)*I)/2))))) + 4)* (1/Decimal(str(np.sinh(np.float128((Decimal(np.pi)*I)/2))))**2) - 4* Decimal((mu**2)*(I**2 - 1)*(I**4)*(Decimal(str(np.sinh(np.float128(Decimal(np.pi)/2))))**2)*(1/Decimal(str(np.sinh(np.float128((Decimal(np.pi)*I)/2))))**2) + 1).sqrt() + 4))
        den = (2*((I**2 - 1)**2)*Decimal((mu**2)*(I**2 - 1)*(I**4)*(Decimal(str(np.sinh(np.float128(Decimal(np.pi)/2))))**2)*(1/Decimal(str(np.sinh(np.float128((Decimal(np.pi)*I)/2))))**2) + 1).sqrt())
        return (num/den)
    else:
        return 0

def abs_der_fun(I,gamma):
    if I!= 0 and I!=1:
        return (np.abs(der_func_f(I,gamma)/Decimal(str((np.sqrt(np.float128(1-function_f(I,gamma)**2)))))))
    elif I== 1:
        return np.abs(np.sinh(np.pi/2)*((2*np.sinh(np.pi/2)**2)*(np.tan(gamma*np.pi/2)**2) + (-np.pi - 4)*np.exp(-np.pi) + 8 + (np.pi - 4)*np.exp(np.pi))*np.pi(tan(gamma*np.pi/2))/(np.sqrt(-np.sinh(np.pi/2)**2*(np.tan(gamma*np.pi/2)**2) + np.exp(np.i) - 2 + np.exp(-np.pi))*(-np.exp(np.pi/2) + exp(-np.pi/2))**2))
    else:
        return 2*np.tan(gamma*np.pi/2)*np.sqrt(np.cosh(np.pi/2)*+2 - 1)/np.pi


a= max_value(abs_der_fun, I_init, gamma)[0]
# I_values = domain_I(I_init, gamma, 0.000001)
# for I in I_values:
#     print("der_f = ", der_func_f(I, gamma), "  abs_der_fun= ", abs_der_fun(I,gamma) , "  I= ", I,"\n")

b =2*Decimal(str(np.tan(np.float128(gamma*Decimal(np.pi)/2))))* Decimal(str(np.sinh(np.float128(Decimal(np.pi/2)))))/Decimal(np.pi)

print(a-b)
# Mu = Decimal(np.tan(np.pi*gamma/2))
# print(Mu)
# #print(np.tan(Decimal(np.pi)*Decimal(gamma)/Decimal(2)))
# print(Decimal(2)*Mu*Decimal(np.sinh(np.pi/2))/Decimal(np.pi))
# conclusion: function f is upper bounded by 1. This implies that \tau is upper
# bounded by pi.

# to find the upper bound of the derivative of tau is necessary to find an
#for its formulat at I = 0 and I= 1.
