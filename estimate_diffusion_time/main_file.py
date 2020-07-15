# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
from estimate_diffusion_time import *
from domain_highways import *
from alpha_function import *
from simpson_rule import *
from integrand import *
import pandas as pd

h = 0.01
interval_eps = np.linspace(0.004 , 0.01 , 4)
interval_angle = np.linspace(0.01 , np.pi/4 - 0.01 , 50)#np.array([0.01 , np.pi/8 , np.pi/4 , 3*np.pi/8 , np.pi/2 - 0.01])
data_time= []

for angle in interval_angle:
    print 'mu = ', 1/np.tan(angle)
    interval_I = def_interval_I(angle)
    print 'domain', interval_I
    #A = alpha_function(I_max)
    for eps in interval_eps:
        list = [angle , eps]
        total_time = 0
        for I in interval_I:
            J= 0.01
            if I != 0:
                A = alpha_function(I)
                T_s_value = simpson_sum( integrand , h, angle , I ,J )
                #print 'T_s = ', T_s_value
                time_of_diffusion = diffusion_time(T_s_value , angle , A, eps)
                total_time = total_time + time_of_diffusion
                J = I
            else:
                time_of_diffusion = 0
            #print diffusion_time
            list.append(time_of_diffusion)
        data_time.append(list)


data_time = np.array(data_time)
#print data_time
#print data_time[:,0]
