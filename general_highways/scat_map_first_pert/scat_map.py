# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
import matplotlib.pyplot as plt
from parametros import *
from functions import *
from find_tau import *
from red_poin_fun import *
from conf_plot import *
#==============================================================================
def SM( domain_I , domain_theta ,step_1 ):
    fig = plt.figure(facecolor= 'white')
    poin_fun_theta , poin_fun_I = [] , []
    dotI_theta , dotI_I , dotTheta_theta , dotTheta_I = [] , [] , [ ] , []
    for I in domain_I:
        print 'I = ' , I
        for theta in domain_theta:
            value_of_tau = assign_tau( I , theta )
            poin_fun_theta.append(red_poi_fun(I, theta , value_of_tau) )
            dotI_theta.append( dot_I( I , theta, value_of_tau) )
        poin_fun_I.append( poin_fun_theta)
        poin_fun_theta = []
        dotI_I.append( dotI_theta  )
        dotI_theta = []
    print len(poin_fun_I)
    SM_table = np.array(poin_fun_I)
    print len(SM_table)
    behavior_I =  np.array(dotI_I)
    #behavior_theta =  dotTheta_I
#
    def estr( x , y ):
        estr = SM_table[ y , x ]
        return estr
    #
    def beh_I( x , y ):
        beh = behavior_I[ y , x ]
        return beh
    #
    y = np.arange(0 , len( domain_I ) , 1 )
    x = np.arange( 0 , len( domain_theta ) , 1 )
    x, y = np.meshgrid(x,y)
    z = estr( x , y )
    plt.contour( x , y , z , 20 , colors = 'blue' , linestyles = 'solid' )
    z_1 = beh_I( x , y )
    plt.contourf( x , y , z_1 , [-1000 , -0.000001 , 0.000001 , 1000]  , colors = ( 'red' , 'white' , 'green' ) , alpha = 0.3)
    plt.contour( x , y , z_1 , levels = [ 0 ] , colors = 'black' , linestyles = 'solid' )
    config_graph( domain_I , domain_theta ,step_1)
    #name_file =  str(par.mu)+ '.svg'
    #plt.savefig(name_file)
