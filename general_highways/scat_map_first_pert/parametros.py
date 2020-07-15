# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
#===============================================================================
# parametros

a_1 = 0.6 # mu
step_1 = 0.1
domain_theta = np.linspace( 0 ,2 * np.pi , 2 * np.pi / step_1)
domain_I = [x for x in np.linspace(-5 , 5, 10/step_1) if x != 0. ]
