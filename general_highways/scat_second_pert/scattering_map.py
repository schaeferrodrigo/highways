# -*- coding: utf-8 -*-
#==============================================================================
import numpy as np
import parameeters as par
import alpha as al
import tau_star as tau
import matplotlib.pyplot as plt
import functions as fun
import newton_method_vv as new
import bifurcation as bif

step_1 = 0.1
domain_theta = np.linspace( 0 ,2 * np.pi , 2 * np.pi / step_1)
domain_I = [x for x in np.linspace(-3 , 3, 6/step_1) if x != 1 ]

def assign_tau( I, theta ):
    value_of_tau_pos = tau.tau( I , theta , 'pos' )
    value_of_tau_neg = tau.tau( I , theta , 'neg' )
    if np.minimum( np.abs(value_of_tau_neg) , np.abs(value_of_tau_pos) ) == np.abs(value_of_tau_neg):
        val_of_tau = value_of_tau_neg
    else:
        val_of_tau = value_of_tau_pos
    return val_of_tau

def config_graph( domain_I , domain_theta , step):
    plt.axis([0,len(domain_theta),0,len(domain_I)])
    plt.xlabel(r'$\theta$', fontsize = 30)
    plt.ylabel('$I$' , fontsize = 30 )
    plt.yticks(np.arange(0,len(domain_I),1/step),('-3','-2','-1','0','1','2','3'),fontsize = 20 )
    plt.xticks(np.arange(0,len(domain_theta),(np.pi/2)/step),('0','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'), fontsize = 20 )

def SM( domain_I , domain_theta ,step_1 ):
    fig = plt.figure(facecolor= 'white')
    bifurcation = bif.singularity_in_I( domain_I , par.mu)
    for I in bifurcation:
        plt.axhline( (I + 3)/step_1 , linestyle = 'dashed' , color = 'black' , linewidth = '1.5' )
    poin_fun_theta , poin_fun_I = [] , []
    dotI_theta , dotI_I , dotTheta_theta , dotTheta_I = [] , [] ,[] ,[]
    tang_min , tang_max , tang_I = [] , [] , []
    for I in domain_I:
        print 'I = ' , I
        for theta in domain_theta:
            value_of_tau = assign_tau( I , theta )
            poin_fun_theta.append( fun.red_poi_fun(I, theta , value_of_tau) )
            dotI_theta.append( fun.dot_I( I , theta, value_of_tau) )
            dotTheta_theta.append( fun.dot_theta( I , theta , value_of_tau ))
        poin_fun_I.append( poin_fun_theta )
        poin_fun_theta = []
        dotI_I.append( dotI_theta )
        dotTheta_I.append( dotTheta_theta)
        dotI_theta , dotTheta_theta = [] , []
    SM_table = np.array(poin_fun_I)
    behavior_I = np.array( dotI_I)
    behavior_theta = np.array( dotTheta_I)
#
    def estr( x , y ):
        estr = SM_table[ y , x ]
        return estr
    #
    def beh_I( x , y ):
        beh = behavior_I[ y , x ]
        return beh
    #
    def beh_theta( x , y ):
        beh = behavior_theta[ y , x ]
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
    z_2 = beh_theta( x , y )
    plt.contour( x , y , z_2 , levels = [ 0 ] , colors = 'red' , linestyles = 'solid' )
    plt.plot( tang_max , tang_I , '.' , color = 'green')
    plt.plot( tang_min , tang_I , '.' , color = 'green')
    config_graph( domain_I , domain_theta ,step_1)
    #fig.text(0.06 , 0.3 , round(par.mu , 2) , ha = 'left', fontsize = 15)
    #fig.text(0.01 , 0.3 , r'$\mu=$', ha = 'left' , fontsize = 15)
    name_file =  str(par.mu)+ '.png'
    plt.savefig(name_file)

SM(domain_I, domain_theta , step_1)
plt.show()
