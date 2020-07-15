# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
import matplotlib.pyplot as plt
#==============================================================================
def config_graph( domain_I , domain_theta , step):
    plt.axis([0,len(domain_theta),0,len(domain_I)])
    plt.xlabel(r'$\theta$', fontsize = 30)
    plt.ylabel('$I$' , fontsize = 30 )
    plt.yticks(np.arange(0,len(domain_I),1/step),('-5','-4','-3','-2','-1','0','1','2','3','4','5'),fontsize = 20 )
    plt.xticks(np.arange(0,len(domain_theta),(np.pi/2)/step),('0','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'), fontsize = 20 )
