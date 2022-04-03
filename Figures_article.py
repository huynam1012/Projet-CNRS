# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 22:17:04 2019

@author: monniello
"""

import numpy as np
import matplotlib.pyplot as plt
from DataArticle import *
import Func_utils as FUN


####################################################
############      FIGURE 2         #################
####################################################

fig2=plt.figure()
ax21 = fig2.add_subplot(211)
ax21.plot(EnJ10,dataJ10[9],'black')
ax21.plot(EnJ10,dataJ10[15],'blue')
ax21.plot(EnJ10,dataJ10[17],'red')
ax21.set_xlim(1.5,3)
ax21.set_ylim(-0.2,0.3)


ax22 = fig2.add_subplot(212)
ax22.plot(EnK10B,dataK10B[9],'black')
ax22.plot(EnK10B,dataK10B[15],'blue')
ax22.plot(EnK10B,dataK10B[17],'red')
ax22.set_xlim(1.5,3)
ax22.set_ylim(-0.2,0.3)
fig2.show()


####################################################
##########      FIGURE 3 a,b,c,d     ###############
####################################################

fig4a=plt.figure()
En = np.arange(1.5,3,0.01)
ax41 = fig4a.add_subplot(2,2,1)
ax41.plot(En,np.real(FUN.RE(1240/En,100)),'black') #The abscisse is the Energy but the function RE takes the wavelength
ax41.plot(En,np.imag(FUN.RE(1240/En,100)),'blue')
ax41.plot(En,FUN.R01(1240/En),'m--')
ax41.plot(En,FUN.R02(1240/En),'g--')
ax41.set_xlim(1.5,3)
ax41.set_ylim(-0.8,0.8)
ax41.set_xlabel('Energy (eV)')

ax42 = fig4a.add_subplot(2,2,3)
ax42.plot(En,np.abs(FUN.RE(1240/En,100))**2,'black')
ax42.plot(En,np.abs(FUN.RE(1240/En,300))**2,'blue')
ax42.plot(En,np.abs(FUN.RE(1240/En,500))**2,'red')
ax42.plot(En,FUN.R01(1240/En)**2,'m--')
ax42.plot(En,FUN.R02(1240/En)**2,'g--')
ax42.set_xlim(1.5,3)
ax42.set_ylim(-0,0.5)
ax42.set_xlabel('Energy (eV)')

ax43 = fig4a.add_subplot(2,2,2)
ax43.plot(En,np.real(FUN.RA(1240/En,100)),'black') #The abscisse is the Energy but the function RE takes the wavelength
ax43.plot(En,np.imag(FUN.RA(1240/En,100)),'blue')
ax43.plot(En,FUN.R01(1240/En),'m--')
ax43.plot(En,FUN.R02(1240/En),'g--')
ax43.set_xlim(1.5,3)
ax43.set_ylim(-1,2)
ax43.set_xlabel('Energy (eV)')

ax44 = fig4a.add_subplot(2,2,4)
ax44.plot(En,np.abs(FUN.RA(1240/En,100))**2,'black')
ax44.plot(En,np.abs(FUN.RA(1240/En,300))**2,'blue')
ax44.plot(En,np.abs(FUN.RA(1240/En,500))**2,'red')
ax44.plot(En,(1+FUN.R01(1240/En))**2,'m--')
ax44.plot(En,(1+FUN.R02(1240/En))**2,'g--')
ax44.set_xlim(1.5,3)
ax44.set_ylim(-0,3)
ax44.set_xlabel('Energy (eV)')

####################################################
##########      FIGURE 4 e,f         ###############
####################################################

def Chi(lambdas,lambda0,Gamma,A,A1,A2):
    j=complex(0,1)
    non_res = A1 + j*A2
    res=FUN.imagLorentz(lambda0,Gamma,A) #res is an complexe lorentzian function. Here we create an "object" containing the function
    return non_res + res.function(lambdas,*res.param) #res.function call the function itself

def contrast(lambdas,delta,delta0,lambda0,Gamma,A,A1,A2,phi,d,f,fy):
    """returns a table with the value of the contrast for each lambda
    """
    chi=Chi(lambdas,lambda0,Gamma,A,A1,A2)
    Re=FUN.RE(lambdas,d);Ra=FUN.RA(lambdas,d)
    cont=[]
    for c,re,ra in zip(chi,Re,Ra):
        function=FUN.Contrast(delta0,c.real,c.imag,f)
        function.fy=fy
        function.RE=re
        function.RA=ra
        function.d=d
        cont.append(function.function(delta,*function.param))
    return np.asarray(cont)

def contrastSiO2(lambdas,delta,delta0,lambda0,Gamma,A,A1,A2,phi,f,fy):
    chi=Chi(lambdas,lambda0,Gamma,A,A1,A2)
    Re=FUN.R01(lambdas);Ra=1+FUN.R01(lambdas)
    cont=[]
    for c,re,ra in zip(chi,Re,Ra):
        function=FUN.Contrast(delta0,c.real,c.imag,f)
        function.fy=fy
        function.RE=re
        function.RA=ra
        cont.append(function.function(delta,*function.param))
    return np.asarray(cont)

def contrastSi(lambdas,delta,delta0,lambda0,Gamma,A,A1,A2,phi,f,fy):
    chi=Chi(lambdas,lambda0,Gamma,A,A1,A2)
    Re=FUN.R02(lambdas);Ra=1+FUN.R02(lambdas)
    cont=[]
    for c,re,ra in zip(chi,Re,Ra):
        function=FUN.Contrast(delta0,c.real,c.imag,f)
        function.fy=fy
        function.RE=re
        function.RA=ra
        cont.append(function.function(delta,*function.param))
    return np.asarray(cont)

fig4b = plt.figure()
ax45 = fig4b.add_subplot(3,2,1)
ax45.plot(En,contrast(1240/En, 0.01, 0, 450, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 80, 0, 0),'blue')
ax45.plot(En,contrast(1240/En, 0.01, 0, 550, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 80, 0, 0),'green')
ax45.plot(En,contrast(1240/En, 0.01, 0, 650, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 80, 0, 0),'orange')
ax45.plot(En,contrast(1240/En, 0.01, 0, 750, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 80, 0, 0),'red')
ax45.set_ylim(-0.005,0.04)

ax46 = fig4b.add_subplot(3,2,3)
ax46.plot(En,contrast(1240/En, 0.01, 0, 450, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 100, 0, 0),'blue')
ax46.plot(En,contrast(1240/En, 0.01, 0, 550, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 100, 0, 0),'green')
ax46.plot(En,contrast(1240/En, 0.01, 0, 650, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 100, 0, 0),'orange')
ax46.plot(En,contrast(1240/En, 0.01, 0, 750, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 100, 0, 0),'red')
ax46.set_ylim(-0.005,0.04)

ax47 = fig4b.add_subplot(3,2,5)
ax47.plot(En,contrast(1240/En, 0.01, 0, 450, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 120, 0, 0),'blue')
ax47.plot(En,contrast(1240/En, 0.01, 0, 550, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 120, 0, 0),'green')
ax47.plot(En,contrast(1240/En, 0.01, 0, 650, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 120, 0, 0),'orange')
ax47.plot(En,contrast(1240/En, 0.01, 0, 750, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 120, 0, 0),'red')
ax47.set_ylim(-0.005,0.04)

ax48 = fig4b.add_subplot(3,2,2)
ax48.plot(En,5*contrastSi(1240/En, 0.01, 0, 450, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 0, 0),'blue')
ax48.plot(En,5*contrastSi(1240/En, 0.01, 0, 550, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 0, 0),'green')
ax48.plot(En,5*contrastSi(1240/En, 0.01, 0, 650, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 0, 0),'orange')
ax48.plot(En,5*contrastSi(1240/En, 0.01, 0, 750, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 0, 0),'red')
ax48.plot(En,contrastSiO2(1240/En, 0.01, 0, 450, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 0, 0),'blue')
ax48.plot(En,contrastSiO2(1240/En, 0.01, 0, 550, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 0, 0),'green')
ax48.plot(En,contrastSiO2(1240/En, 0.01, 0, 650, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 0, 0),'orange')
ax48.plot(En,contrastSiO2(1240/En, 0.01, 0, 750, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 0, 0),'red')
ax48.set_ylim(-0.04,0.005)

ax49 = fig4b.add_subplot(3,2,4)
ax49.plot(En,contrast(1240/En, 0.01, 0, 450, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 300, 0, 0),'blue')
ax49.plot(En,contrast(1240/En, 0.01, 0, 550, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 300, 0, 0),'green')
ax49.plot(En,contrast(1240/En, 0.01, 0, 650, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 300, 0, 0),'orange')
ax49.plot(En,contrast(1240/En, 0.01, 0, 750, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 300, 0, 0),'red')
ax49.set_ylim(-0.005,0.04)

ax410 = fig4b.add_subplot(3,2,6)
ax410.plot(En,contrast(1240/En, 0.01, 0, 450, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 500, 0, 0),'blue')
ax410.plot(En,contrast(1240/En, 0.01, 0, 550, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 500, 0, 0),'green')
ax410.plot(En,contrast(1240/En, 0.01, 0, 650, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 500, 0, 0),'orange')
ax410.plot(En,contrast(1240/En, 0.01, 0, 750, 0.00004, 0.000000001, 0.00001, 0.00001, 1/2, 500, 0, 0),'red')
ax410.set_ylim(-0.005,0.04)













