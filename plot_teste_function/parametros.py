# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
#===============================================================================
# parametros

gamma = 0.4 # gamma domain (0,0,5)
a_1 = np.sin(gamma * np.pi/2)
a_2 = np.cos(gamma * np.pi/2)
mu  = a_1/a_2
step_1 = 0.1
domain_theta = np.linspace( 0 ,2 * np.pi , 10000)
r = 0.01 # domain of r [0,1]
I = 2.
