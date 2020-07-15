# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
import matplotlib.pyplot as plt

def alpha_function(I):
    if I != 0 :
        alpha = np.sinh(np.pi/2)* (I**2)/np.sinh(np.pi* I/2)
    else:
        alpha = 0
    return alpha

# alpha_values = []
# interval_I= np.linspace(0.01 , 2 , 1.99/0.01)
# for I in interval_I:
#     alpha_value = alpha_function(I)
#     alpha_values.append(alpha_value)
#
# plt.plot(interval_I , alpha_values)
# plt.show()
