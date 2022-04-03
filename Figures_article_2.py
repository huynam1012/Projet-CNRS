# -*- coding: utf-8 -*-
"""
Created on Sun Jan 20 12:52:03 2019

@author: monniello
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from DataArticle import *
import Func_utils as FUN


####################################################
##########      FIGURE 5             ###############
####################################################

#First we fit the data to obtain the parameters

function1 = [FUN.Polynom(2,[10000,1,0.5])]
function2 = [FUN.Polynom(2,[1000,1,0.])]
function3 = [FUN.Polynom(2,[1000,1,0])]

Fit1 = FUN.plot_fit(AnglesPolJ10,function1,[polarJ10[10]/3000],[0])
Fit2 = FUN.plot_fit(AnglesPolJ10,function2,[polarJ10[500]/3000],[0])
Fit3 = FUN.plot_fit(AnglesPolJ10,function3,[polarJ10[1000]/3000],[0])


delta0=[]
f=[]
function=[FUN.Polynom2(10000,0,0.5)]
FIT=FUN.plot_fit(AnglesPolJ10,function,polarJ10,shift=False)
delta0=FIT[1].T[1]*1000
f=FIT[1].T[2]*10000


#Now we plot the data

fig = plt.figure(figsize=(8,3))
gs=GridSpec(2,2) #Allows to have figures on a specific display

ax1=fig.add_subplot(gs[:,0]); #first figure is on 2 rows and 1 column
ax1.plot(EnJ10,polarJ10.T[9]/3000,'black')
ax1.plot(EnJ10,polarJ10.T[14]/3000,'blue')
ax1.plot(EnJ10,polarJ10.T[21]/3000,'green')
ax1.plot(EnJ10,polarJ10.T[27]/3000,'orange')
ax1.plot(EnJ10,polarJ10.T[32]/3000,'red')
ax1.set_xlim(1.5,3)
ax1.set_ylim(0,10)
ax1.set_xlabel('Energy (eV)')
ax1.set_ylabel('Intensity (arb. unit)' )

ax2 = fig.add_subplot(gs[1,1])
ax2.plot(EnJ10,delta0)
ax2.set_xlim(1.5,3)
ax2.set_ylim(-3,3)

ax3 = fig.add_subplot(gs[0,1])
ax3.plot(EnJ10,f)
ax3.set_xlim(1.5,3)
ax3.set_ylim(0,0.8)
fig.tight_layout()

####################################################
##########         Inset  Figure 5        ##########
####################################################

x=np.arange(-0.04,0.04,0.002)

subpos = [0.07,0.6,0.4,0.3] #Positions of the inset are defined in relative : 0.1 is 10% of the size
subax1 = FUN.add_subplot_axes(ax1,subpos)
subax1.plot(AnglesPolJ10,polarJ10[10]/3000,'blue',linestyle = 'none', marker='o', markersize=3, fillstyle='none')
function1[0].draw(x,'blue',axes=subax1)
subax1.plot(AnglesPolJ10,polarJ10[500]/3000,'green',linestyle = 'none',marker='o', markersize=3, fillstyle='none')
function2[0].draw(x,'green',axes=subax1)
subax1.plot(AnglesPolJ10,polarJ10[1000]/3000,'red',linestyle = 'none',marker='o', markersize=3, fillstyle='none')
function3[0].draw(x,'red',axes=subax1)
subax1.set_xlim(-0.04,0.04)
subax1.set_ylim(-0.3,3)


####################################################
##########      FIGURE 6             ###############
####################################################


####################################################
##########      Tube A  K10          ###############
####################################################
dataK=lambdasK10,dataK10,AnglesK10

FIT_K10_ideal_function=FUN.fit_fig6(*dataK,FUN.Contrast_ideal(0.0001,.00001))
FIT_K10_NC_function=FUN.fit_fig6(*dataK,FUN.Contrast_NCoherent(0.0001,0.0001,0.01))
FIT_K10_function=FUN.fit_fig6(*dataK,FUN.Contrast(-0.001,0.000005,0.00005,0.01))

####  Extracttion of delta0 and f from the data of thius tube
    
delta0K10=[]
fK10=[]
function=[FUN.Polynom2(10000,0,0.5)]
FIT=FUN.plot_fit(AnglesPolK10,function,polarK10,shift=False)
delta0K10=FIT[1].T[1]*1000
fK10=FIT[1].T[2]*10000    

####################################################
##########      Tube B  J10          ###############
####################################################
dataJ=lambdasJ10,dataJ10,AnglesJ10

FIT_J10_ideal_function=FUN.fit_fig6(*dataJ,FUN.Contrast_ideal(0.0001,.00001))
FIT_J10_NC_function=FUN.fit_fig6(*dataJ,FUN.Contrast_NCoherent(0.00001,0.00001,0.001))
FIT_J10_function=FUN.fit_fig6(*dataJ,FUN.Contrast(-0.001,0.000005,0.00005,0.01))

####  Extracttion of delta0 and f from the data of this tube
    
delta0J10=[]
fJ10=[]
function=[FUN.Polynom2(10000,0,0.5)]
FIT=FUN.plot_fit(AnglesPolJ10,function,polarJ10,shift=False)
delta0J10=FIT[1].T[1]*1000
fJ10=FIT[1].T[2]*10000    

################################################
#############   Figure                   #######
################################################

fig6 = plt.figure(figsize=(8,6))
gs=GridSpec(3,2,height_ratios=(3,1,1))
n=630

ax1 = fig6.add_subplot(gs[0,0])
ax1.plot(dataK[2],dataK[1].T[n],linestyle = 'none', marker='o',fillstyle='none')
FIT_K10_ideal_function[n][0].draw(np.arange(AnglesK10[0],AnglesK10[-1],-0.0001),'black',axes=ax1)
FIT_K10_NC_function[n][0].draw(np.arange(AnglesK10[0],AnglesK10[-1],-0.0001),'orange',axes=ax1)
FIT_K10_function[n][0].draw(np.arange(AnglesK10[0],AnglesK10[-1],-0.0001),'red',axes=ax1)
FIT_K10_function[n][0].draw_scat(np.arange(AnglesK10[0],AnglesK10[-1],-0.0001),'green',axes=ax1)
ax1.set_xlim(-0.1,0.1)
ax1.set_ylim(-0.1,0.2)

#figure delta0
ax2 = fig6.add_subplot(gs[2,0])
ax2.plot(EnK10B,delta0K10,'blue')
ax2.plot(EnK10B,[1000*function[0].param[0] for function in FIT_K10_function[1:]],'red')
ax2.set_xlim(1.5,3)
ax2.set_ylim(-4,4)

#figure f
ax3 = fig6.add_subplot(gs[1,0])
ax3.plot(EnK10B,10*np.sqrt(abs(fK10)),'blue')
ax3.plot(EnK10B,[1000*abs(function[0].param[2]) for function in FIT_K10_NC_function[1:]],'orange')
ax3.plot(EnK10B,[1000*function[0].param[3] for function in FIT_K10_function[1:]],'red')
ax3.set_xlim(1.5,3)
ax3.set_ylim(-0,15)

n=630
ax4 = fig6.add_subplot(gs[0,1])
ax4.plot(dataJ[2],dataJ[1].T[n],linestyle = 'none', marker='o',fillstyle='none',axes=ax4)
FIT_J10_ideal_function[n][0].draw(np.arange(AnglesJ10[0],AnglesJ10[-1],-0.0001),'black',axes=ax4)
FIT_J10_NC_function[n][0].draw(np.arange(AnglesJ10[0],AnglesJ10[-1],-0.0001),'orange',axes=ax4)
FIT_J10_function[n][0].draw(np.arange(AnglesJ10[0],AnglesJ10[-1],-0.0001),'red',axes=ax4)
FIT_J10_function[n][0].draw_scat(np.arange(AnglesJ10[0],AnglesJ10[-1],-0.0001),'green',axes=ax4)
ax4.set_xlim(-0.1,0.1)
ax4.set_ylim(-0.1,0.2)

#figure delta0
ax5 = fig6.add_subplot(gs[2,1])
ax5.plot(EnJ10,delta0J10,'blue')
ax5.plot(EnJ10,[1000*function[0].param[0] for function in FIT_J10_function[1:]],'red')
ax5.set_xlim(1.5,3)
ax5.set_ylim(-4,4)

#figure f
ax6 = fig6.add_subplot(gs[1,1])
ax6.plot(EnJ10,10*np.sqrt(abs(fJ10)),'blue')
ax6.plot(EnJ10,[1000*abs(function[0].param[2]) for function in FIT_J10_NC_function[1:]],'orange')
ax6.plot(EnJ10,[1000*function[0].param[3] for function in FIT_J10_function[1:]],'red')
ax6.set_xlim(1.5,3)
ax6.set_ylim(-0,15)




######################################################
#########              Figure 7           ############
######################################################

Real_K10 = [2000*function[0].param[1] for function in FIT_K10_function[1:]]
Imag_K10 = [2000*function[0].param[2] for function in FIT_K10_function[1:]]
KK_K10 = FUN.Kramers(1240/lambdasK10[::-1],Imag_K10[::-1],0.001)

Real_J10 = [2000*function[0].param[1] for function in FIT_J10_function[1:]]
Imag_J10 = [2000*function[0].param[2] for function in FIT_J10_function[1:]]
KK_J10 = FUN.Kramers(1240/lambdasJ10[::-1],Imag_J10[::-1],0.001)


plt.figure()
plt.plot(1240/dataK[0],Imag_K10)
plt.plot(1240/dataK[0],Real_K10)
plt.plot(np.arange(1.5,3,0.001),KK_K10)
plt.ylim(-1,1)
plt.xlim(1.5,3)


plt.figure()
plt.plot(1240/dataJ[0],Imag_J10)
plt.plot(1240/dataJ[0],Real_J10)
plt.plot(np.arange(1.5,3,0.001),KK_J10)
plt.ylim(-1,1)
plt.xlim(1.5,3)






