# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
import bigfloat as bf
#===============================================================================
# parametros

gamma= 0.4 # gamma domain (0,0,5)
a_1 = bf.sin(gamma * np.pi/2)
a_2 = bf.cos(gamma * np.pi/2)
print( "mu=",a_1/a_2)
step_1 = 0.1
domain_theta = np.linspace( 0 ,2 * np.pi , 70)
r = 0.01 # domain of r [0,1]
domain_I = [x for x in np.linspace(-5 , 5, 100) if x != 0. ]
