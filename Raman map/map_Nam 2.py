# -*- coding: utf-8 -*-

from colorama import Fore, Back, Style 
import matplotlib.pyplot as plt
import numpy as np
from math import *
import peakutils
from lmfit.models import GaussianModel, LinearModel, LorentzianModel
from lmfit.models import GaussianModel, LinearModel, LorentzianModel, BreitWignerModel
from progressbar import ProgressBar
pbar = ProgressBar()
from matplotlib import pyplot
import peakutils
from scipy import stats
import os 
import pylab

RBM_list=[]
G_list=[]
filename = 'NT3_power1500microw_10s.txt' 
nb_pics=2
DATA=np.loadtxt(filename)
line=np.shape(DATA)[0]
G_band_pos=[]
diameter_list=[]
rapport_G_D_LIST=[]
cell_counter=1
line_counter=1
for a in range(1,line):
    if DATA[a-1,2]<DATA[a,2]:
        cell_counter +=1
    if DATA[a-1,1] != DATA[a,1]:
        line_counter +=1
C=int(cell_counter/line_counter)
L=int(line_counter)
MAP=np.zeros((L,C))
NT_counter=0
cible=1592
element=int(line/cell_counter)
for a in pbar(range(L)):
    for b in range(C):
        adresse_in=a*C*element+b*element
        adresse_out=a*C*element+b*element+element
        X=DATA[adresse_in:adresse_out,2]
        Y=DATA[adresse_in:adresse_out,3]
        MAP[a,b]=0
        G_OUT=int(np.argmin(np.abs(X-cible+100)))
        G_IN=int(np.argmin(np.abs(X-cible-100)))                
        NT_counter+=1

        x=X[G_IN:G_OUT]
        y=Y[G_IN:G_OUT]
        

        #name=str(a)+str('-')+str(b)
        #pylab.savefig(name+str('.png'))
        
        moyenne_locale=np.mean(y)
        if(Y[int(np.argmin(np.abs(X-cible)))])>3*moyenne_locale:
        
            #background = LinearModel(prefix='bck_')
            #pars = background.guess(y, x=x)
                    
            #M1 = LorentzianModel(prefix='M1_')
            #pars.update(M1.make_params(center=1590))
                    
            #M2 = LorentzianModel(prefix='M2_')
            #pars.update(M2.make_params(center=1550))
            #pars['M2_center'].set(min=1500,max=1600)
                    
            #mod=background+M1+M2
            #init = mod.eval(pars, x=x)
            #out = mod.fit(y, pars, x=x)
            #comps = out.eval_components(x=x)
    
    
                    
            #model=comps['bck_']+comps['M1_']+comps['M2_']
            #error=np.std(np.abs(model-y)/np.mean(model))
    
          
            #plt.plot(x, out.best_fit, 'r-',x,y)
            plt.plot(x,y)
            
            #print(error)
            plt.show()
            #sauvegarde
            MAP[a,b]=Y[int(np.argmin(np.abs(X-cible)))]
            
plt.imshow(MAP,cmap="nipy_spectral")       
plt.colorbar()