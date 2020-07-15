# -*- coding: utf-8 -*-
#==============================================================================
import numpy as np
import functions as fun
import matplotlib.pyplot as plt
import parameeters as par
import alpha as al

def eixo_theta( phi ):
    eixo = phi
    return phi

def plot_crista( I , domain_angle):
    ''' plot das cristas'''
    fig = plt.figure()
    if np.abs( al.alpha( I ) * par.mu ) < 1:
        plt.plot( domain_angle , fun.xi_hor( I , domain_angle ) , '-' , color =  'red'  )
        plt.plot( domain_angle , -fun.xi_hor( I , domain_angle) + np.pi , '-' , color = 'blue' )
        plt.plot( domain_angle , fun.xi_hor( I , domain_angle) + 2*np.pi , '-' , color = 'red' )
        plt.axis([0,2*np.pi,-np.pi/2 , 5*np.pi/2])
        plt.yticks(np.arange( -np.pi/2 , 5*np.pi/2 , np.pi/2 ),('$-\pi/2$' , '0', '$\pi/2$','$\pi$','$3\pi/2$' , '$2\pi$' , '$5\pi/2$' ) , fontsize = 25 )
        plt.xticks(np.arange(0,2*np.pi,np.pi/2),('0','$\pi/2$','$\pi$','$3\pi/2$'), fontsize = 25)
    else:
        plt.plot( fun.xi_vert(I , domain_angle) , domain_angle , '-' , color = 'red')
        plt.plot( -fun.xi_vert(I , domain_angle) + np.pi , domain_angle , '-' , color = 'blue')
        plt.plot( fun.xi_vert(I , domain_angle) + 2*np.pi, domain_angle , '-' , color = 'red')
        plt.axis([-np.pi/2 , 5*np.pi/2 , 0 , 2*np.pi ])
        plt.yticks(np.arange(0,2*np.pi,np.pi/2),('0','$\pi/2$','$\pi$','$3\pi/2$'), fontsize = 25)
        plt.xticks(np.arange( -np.pi/2 , 5*np.pi/2 , np.pi/2 ),('$-\pi/2$' , '0', '$\pi/2$','$\pi$','$3\pi/2$' , '$2\pi$' , '$5\pi/2$' ) , fontsize = 25 )
    plt.plot( domain_angle , eixo_theta( domain_angle ) , '-' , color = 'black')
    plt.xlabel(r'$\varphi$', fontsize = 30)
    plt.ylabel( r'$\sigma$' , fontsize = 30 , rotation = 'horizontal' )
    fig.text(0.06 , 0.05 , round(par.mu , 2) , ha = 'left', fontsize = 15)
    fig.text(0.01 , 0.05 , r'$\mu=$', ha = 'left' , fontsize = 15)
    fig.text(0.06 , 0.01 , I , ha = 'left', fontsize = 15)
    fig.text(0.01 , 0.01 , r'$I$', ha = 'left' , fontsize = 15)
    plt.show()

domain_phi = np.arange(0,2*np.pi,0.01)
plot_crista( 0.9935 , domain_phi)
