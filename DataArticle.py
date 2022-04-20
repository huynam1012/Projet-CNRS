# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 20:57:09 2019

@author: monniello
"""

import numpy as np
import Func_utils as FUN


def n1(lambdas):
    a1=1
    a2=0.6961663/(1-(0.0684043/(0.001*lambdas))**2)
    a3=0.4079426/(1-(0.1162414/(0.001*lambdas))**2)
    a4=0.8974794/(1-(9.896161/(0.001*lambdas))**2)    
    return np.sqrt(a1+a2+a3+a4)

def n2(lambdas):
    return 3.555 + 15*np.exp(-0.001*lambdas/0.16452)+884775.04994*np.exp(-0.001*lambdas/0.02857)

def R01(lambdas):
    return (1-n1(lambdas))/(1+n1(lambdas))

def R12(lambdas):
    return (n1(lambdas)-n2(lambdas))/(n1(lambdas)+n2(lambdas))

def RE(lambdas,d):
    j=complex(0,1)
    num = R01(lambdas) + R12(lambdas)*np.exp(2*1.48*np.pi*j*2*d/lambdas)
    denom = 1 + R01(lambdas)*R12(lambdas)*np.exp(2*1.48*np.pi*j*2*d/lambdas)
    return num/denom

def RA(lambdas,d):
    j=complex(0,1)
    num=(1-R01(lambdas)**2)*np.exp(2*1.48*np.pi*j*2*d/lambdas)
    denom=1 + R01(lambdas)*R12(lambdas)*np.exp(2*1.48*np.pi*j*2*d/lambdas)
    return 1-num/denom

repertoire = 'C:/Users/Nam/Desktop/Article Montpellier/'

#########################################
#########       TUBEK10           #######
#########################################
    
name = 'contrast_smoothed_20170310.dat'
dataK10 = np.genfromtxt(repertoire+name).T

lambdasK10=dataK10[0]
EnK10 = dataK10[1]
dataK10=dataK10[2:]

AnglesK10 = np.asarray([-15,-10,-6,-3,-1.5,-1,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,0,0.2,0.5,2,5,9,14])
AnglesK10 = -AnglesK10*np.pi/180-0.001

polarK10 = np.genfromtxt(repertoire+'polar_20170310.dat')
polarK10=polarK10.T[2:].T 
AnglesPolK10 = np.asarray([-((-2+0.1*(i))*np.pi/180+0.001) for i,_ in enumerate(polarK10[0])])


#########################################
#########       TUBEK10_bis       #######
#########################################
    
name = 'contrast_smoothed_20170321_K10.dat'
dataK10B = np.genfromtxt(repertoire+name).T

lambdasK10B=dataK10B[0]
EnK10B = dataK10B[1]
dataK10B=dataK10B[2:]

polarK10B=np.genfromtxt(repertoire + 'polar20170321_K10.dat')
polarK10B=polarK10B.T[2:].T
AnglesPolK10B = np.asarray([-((-2+0.1*(i))*np.pi/180+0.001) for i,_ in enumerate(polarK10B[0])])

AnglesK10B = np.asarray([-15,-10,-6,-3,-1.5,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.2,0.5,2,5,9,14])
AnglesK10B = -AnglesK10B*np.pi/180-0.001


#########################################
#########       TUBEJ10           #######
#########################################
polarJ10=np.genfromtxt(repertoire + 'polar20170321.dat')
polarJ10=polarJ10.T[2:].T
AnglesPolJ10 = np.asarray([-((-2+0.1*(i))*np.pi/180+0.001) for i,_ in enumerate(polarJ10[0])])


name = 'contrast_smoothed_20170321.dat'
dataJ10 = np.genfromtxt(repertoire+name).T

lambdasJ10=dataJ10[0]
EnJ10 = dataJ10[1]
dataJ10=dataJ10[2:]

AnglesJ10 = np.asarray([-15,-10,-6,-3,-1.5,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.2,0.5,2,5,9,14])
AnglesJ10 = -AnglesJ10*np.pi/180-0.001

