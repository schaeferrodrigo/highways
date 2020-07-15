# -*- coding: utf-8 -*-
#==============================================================================
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
import parameeters as par
import functions as fun

def retas(I , phi):
    reta = ((I-1)/I) * phi
    return reta

def Meln_Pot( I , phi , s):
    '''Potencial de Melnikov'''
    pot = fun.A_1(I) * np.cos( phi ) + fun.A_2(I) * np.cos(phi - s)
    return pot

def plot_cristas_retas( I , phi ):
    ''' plot de cristas e retas'''
    plt.plot( phi , fun.xi_hor( I , phi),  '.' , color = 'red' )
    plt.plot( phi , fun.xi_hor(I , phi )  + 2 * np.pi , '.' , color = 'red' , linewidth = 3 )
    plt.plot( phi , -fun.xi_hor( I ,phi ) + np.pi,'.', color = 'blue' , linewidth = 3 )
    const = 0
    while const < 4*np.pi:
        plt.plot( phi , retas(I,phi) + const, '.' , color = 'black' , markersize = .3 )
        const = const + 0.5

def threed_plot(I , variable):
    '''plot 3D do Potencial de Melnikov
    variable - Escolha de plot em s ou sigma'''
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x = np.arange(0,2*np.pi,0.1)
    phi, s = np.meshgrid(x, x)
    z = Meln_Pot(I, phi , s )
    if variable == 's':
        var = s
        ax.set_ylabel(r'$s$ ',fontsize = 20)
    elif variable == 'sigma':
        var = phi -s
        ax.set_ylabel(r'$\sigma$ ',fontsize = 20)
    cset = ax.contour(phi, var, z, zdir='z', offset = -9, cmap=cm.coolwarm)
    surf = ax.plot_surface(phi, var, z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
    ax.set_zlim(-9,8.5)
    ax.set_xlabel(r'$\varphi$',fontsize = 20)
    ax.zaxis.set_major_locator(LinearLocator(5))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

def twoD_plot( I , variable):
    fig = plt.figure( )
    x = np.arange(0,2*np.pi,0.01)
    if variable == 's':
        phi, s = np.meshgrid(x, x)
        var , text = s ,variable
    elif variable == 'sigma':
        phi,sigma = np.meshgrid(x, x)
        var, text , s = sigma , r'$\sigma$' , phi - sigma
        plot_cristas_retas( I , phi )
    z = Meln_Pot(I, phi , s )
    plt.contour( phi , var , z , 30 )
    plt.xlabel(r'$\varphi$', fontsize = 30)
    plt.ylabel( text , fontsize = 30 , rotation = 'horizontal' )
    plt.axis([0,2*np.pi,0,2*np.pi])
    plt.yticks(np.arange(np.pi/2,2*np.pi,np.pi/2),('$\pi/2$','$\pi$','$3\pi/2$'),fontsize = 30)
    plt.xticks(np.arange(0,2*np.pi,np.pi/2),('0','$\pi/2$','$\pi$','$3\pi/2$'), fontsize = 30)
    plt.show()

twoD_plot( 0.5 , 'sigma')
