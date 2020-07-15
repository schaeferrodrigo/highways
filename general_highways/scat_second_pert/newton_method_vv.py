# -*- coding: utf-8 -*-
#==============================================================================
import numpy as np
import functions as fun
import alpha as al
import parameeters as par
import matplotlib.pyplot as plt
import math
import sys

def F( I , theta , tau ):
    ''' Critério para tangência.
    F1 - equação da crista
    F2 - d^2L^*/dtau^2 = 0 - interseccao nao transversal entre as variedades invariantes'''
    F1 = I * fun.A_1( I ) * np.sin( theta - I * tau ) + ( I - 1 ) * fun.A_2( I ) * np.sin( theta - ( I - 1 ) * tau )
    F2 = - ( I**2 ) * fun.A_1( I ) * np.cos( theta - I * tau ) - ( (I - 1 )**2)* fun.A_2( I ) * np.cos( theta - (I - 1) * tau)
    return np.array( [ [ F1 ] , [ F2 ] ])

def DF( I , theta , tau ):
    ''' Matriz Jacobiana da funçao acima'''
    DF1_theta  = I * fun.A_1( I ) * np.cos( theta - I * tau ) + ( I - 1 )*fun.A_2(I) * np.cos( theta - (I - 1) * tau)
    DF1_tau =  -(I**2) * fun.A_1( I ) * np.cos( theta - I * tau ) - ((I - 1)**2)* fun.A_2( I ) * np.cos( theta - (I - 1) * tau)
    DF2_theta = (I**2) * fun.A_1(I) * np.sin( theta - I * tau ) + (( I - 1 )**2) * fun.A_2(I) * np.sin( theta - (I -1)*tau )
    DF2_tau = -(I**3) * fun.A_1( I ) * np.sin( theta - I * tau ) - ((I- 1)**3) * fun.A_2(I) * np.sin(theta - (I - 1 ) * tau)
    return np.array([[DF1_theta , DF1_tau] , [DF2_theta , DF2_tau]])

def newton_method( I , initial_value, psi , tol = 1e-6):
    '''método de newton para encontrar theta onde ocorre as tangências'''
    theta, tau  = initial_value[ 0 ] , initial_value[ 1 ]
    if np.linalg.norm( F( I , theta, tau ) ) < tol:
        roots = [ theta , tau]
    else:
        while np.linalg.norm(F( I , theta, tau )) > tol:
            vector = np.array([ [ theta ] , [ tau ]])
            newton_step = vector - np.linalg.inv( DF( I , vector[ 0 , 0 ] , vector[ 1 , 0 ] )).dot( F( I , vector[0,0] , vector[ 1 , 0 ] ) )
            theta = ( newton_step[ 0 , 0 ] )%(2 * np.pi)
            tau = newton_step[ 1 , 0 ]
            roots = [ theta , tau]
    return roots

def theta_tangent( I , index ):
    ''' funçao que retorna theta e tau onde corre as tangências'''
    psi = fun.psi( I )[ index ]
    theta_tau = fun.theta( I , psi )
    theta , tau = theta_tau[ 0 ] , theta_tau[ 1 ]
    initial_value = [theta , tau ]
    solution = newton_method( I , initial_value , psi )
    return solution
