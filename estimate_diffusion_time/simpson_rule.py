# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np

def simpson_sum(f , h , angle ,I_f , I_0 = 0.01):
    "Simpson's rule, ref_book = útiles básicos de cálculo numérico"
    net = np.linspace(I_0 , I_f , (I_f - I_0)/h )
    sum = 0
    for i in range(0 , len(net)-1):
        if i%2==1:
            sum = sum + 4 * f(net[i],angle)
        else:
            if i==0 or  i == len(net)-1:
                sum = sum + f(net[i],angle)
            else:
                sum = sum + 2*f(net[i],angle)
    return sum
