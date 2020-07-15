# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
from functions import *
import matplotlib.pyplot as plt
from parametros import *
from red_poin_fun import red_poin_approx
import statistics as st
import bigfloat as bf
#===============================================================================

def function_f(I):
    num = (I**2)*A(-1,a_2)-bf.sqrt(A(-1,a_2)**2 +(I**2-1)*(I**2)*A(I,a_1)**2)
    den = (I**2 - 1) * A(-1,a_2)
    return num/den

def tau_highways(I):
    tau_h = np.sign(I)*bf.acos(function_f(I))
    tau_H = -np.sign(I)*bf.acos(function_f(I))
    return [tau_h , tau_H]

def theta_h(I):
    num,den = A(-1,a_2)*(1-function_f(I)),A(I,a_1)
    main_part = bf.acos(num/den)
    theta = main_part + np.sign(I)*I*bf.acos(function_f(I))
    return theta

def theta_H(I):
    num,den = A(-1,a_2)*(1-function_f(I)), A(I,a_1)
    main_part = bf.acos(num/den)
    theta = -main_part -np.sign(I) *I*bf.acos(function_f(I))
    return (theta)%(2*np.pi)


# def level_curves():
#     highway_h, highway_H = [],[]
#     for I in domain_I:
#         tau_h,tau_H = tau_highways(I)
#         print("tau_h",tau_h)
#         red_Poin_r_h = red_poin_approx(I, theta_h(I),tau_h)
#         red_Poin_r_H = red_poin_approx(I, theta_H(I),tau_H)
#         highway_h.append(red_Poin_r_h)
#         highway_H.append(red_Poin_r_H)
#     m_h,med_h,M_h = min(highway_h),st.median(highway_h),max(highway_h)
#     m_H,med_H,M_H = min(highway_H),st.median(highway_H),max(highway_H)
#     return [m_h,M_h,m_H , M_H]

# construction zone=============================================================
#
# print(level_curves())





# ====================== Test if tau is well defined ===========================
# def theta_h_tau(I):
#     tau_h = tau_highways(I)[0]
#     num,den = A(-1,a_2)*(1-np.cos(tau_h)),A(I,a_1)
#     main_part = np.arccos(num/den)
#     theta = main_part + np.sign(I)*I*np.arccos(np.cos(tau_h))
#     return theta
#
# def theta_H_tau(I):
#     tau_H = tau_highways(I)[1]
#     num,den = A(-1,a_2)*(1-np.cos(tau_H)), A(I,a_1)
#     main_part = np.arccos(num/den)
#     theta = -main_part -np.sign(I) *I*np.arccos(np.cos(tau_H))
#     return (theta)%(2*np.pi)



#======================= Plot of the highways using I and theta defined above
# domain_I = np.array(domain_I)
# #
# plt.plot( theta_H(domain_I),domain_I, "blue")
# plt.plot(theta_h(domain_I),domain_I , "blue")
# plt.plot(theta_H_tau(domain_I), domain_I, "red")
# plt.plot(theta_h_tau(domain_I), domain_I,"red")
# plt.axis([0,2*np.pi,-5,5])
# plt.show()
#
# print( A(-1,a_2)*(1-function_f(domain_I)))
# print(A(domain_I,a_1))
# print(theta_H(domain_I))
