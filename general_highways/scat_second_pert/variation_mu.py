# -*- coding: utf-8 -*-
#==============================================================================
import numpy as np
import parameeters as par
import scattering_map as scat
import matplotlib.pyplot as plt

domain_mu = np.linspace( 0.01 , 0.1 , 0.9/0.009 )
for mu in domain_mu:
    par.mu = mu
    par.a_1 = mu
    step_1 = 0.005
    domain_theta = np.linspace( 0 ,2 * np.pi , 2 * np.pi / step_1)
    domain_I = [x for x in np.linspace(-3 , 3, 6/step_1) if x != 1 ]
    scat.SM(domain_I, domain_theta , step_1)
