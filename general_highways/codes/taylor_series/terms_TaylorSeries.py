# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
from find_tau import *
from red_poin_fun import *
from approx_frp import *
from parametros import *
#===============================================================================
# h = 0.01 # step
# num_terms = 100 # number of terms in the Taylor's series
#
#
# I = 0.1
# theta = theta_H(I)

def taylor_series_value(I, theta, num_terms , h):
    i = 0
    Taylor_series = 0
    coeff_taylor = []
    while i <= num_terms:
        print("series order =", i)
        ith_derivate = 0
        for j in np.arange(0, i+1):
            print("j=" ,j)
            binomial = (np.math.factorial(i)/
                       (np.math.factorial(j)*np.math.factorial(i-j)))
            quotient , r_step  = h**i , (i/2 - j)*h
            tau = assign_tau(I,theta ,r_step)
            RedPoin = red_poi_fun(I, theta , r_step)
            ith_derivate = ((-1)**j) * binomial * RedPoin /quotient + ith_derivate
        ith_Taylor = (r**i) * ith_derivate /np.math.factorial(i)
        Taylor_series = Taylor_series + ith_Taylor
        print(Taylor_series)
        coeff_taylor.append(ith_Taylor)
        i = i+1
    return Taylor_series



def energies( I, theta, num_terms, h):
    taylor_series = taylor_series_value(I,theta,num_terms,h)
    highway = A(-1,a_2)
    tau = assign_tau(I,theta)
    simple_app = red_poin_approx(I,theta,tau)
    return [taylor_series,highway]


#print(energies(0.1,theta_H(0.1), 100 , 0.01))







# i = 0
# Taylor_series = 0
# coeff_taylor = []
# while i <= num_terms:
#     print("series order =", i)
#     ith_derivate = 0
#     for j in np.arange(0, i+1):
#         print("j=" ,j)
#         binomial = (np.math.factorial(i)/
#                    (np.math.factorial(j)*np.math.factorial(i-j)))
#         quotient , r_step  = h**i , (i/2 - j)*h
#         tau = assign_tau(I,theta ,r_step)
#         RedPoin = red_poi_fun(I, theta , r_step)
#         ith_derivate = ((-1)**j) * binomial * RedPoin /quotient + ith_derivate
#     ith_Taylor = (r**i) * ith_derivate /np.math.factorial(i)
#     Taylor_series = Taylor_series + ith_Taylor
#     print(Taylor_series)
#     coeff_taylor.append(ith_Taylor)
#     i = i+1
#
# print("taylor series =",Taylor_series)
#
# print("Highway = ", A(-1,a_2))
#
# tau = assign_tau(I,theta)
# simple_app = red_poin_approx(I,theta,tau)
# print('simple =', simple_app)
