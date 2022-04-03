#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 10:42:54 2018

@author: leonard

Module for fitting data.

contains class functions as Gaussian, Fano, Lorentz, and Polynom.
The purpose of thess classes is to use the function to fit data.
The parameters are stored to serve as seed for the fit function and are to be replaced by the output of the fit after it


"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from functools import partial


# Print iterations progress
# Just to know for some long calculations where we are
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()

def add_subplot_axes(ax,rect,axisbg='w'):
    """function to add a subplot as an insert
    https://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib
    """
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height])#,axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax      

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

def R02(lambdas):
    return (1-n2(lambdas))/(1+n2(lambdas))

def RE(lambdas,d):
    j=complex(0,1)
    num = R01(lambdas) + R12(lambdas)*np.exp(2*1.48*np.pi*j*2*d/lambdas)
    denom = 1 + R01(lambdas)*R12(lambdas)*np.exp(2*1.48*np.pi*j*2*d/lambdas)
    return num/denom

def RA(lambdas,d):
    #In the article R coefficient are taken positive. Here they are negative
    j=complex(0,1)
    num=(1-R01(lambdas)**2)*R12(lambdas)*np.exp(2*1.48*np.pi*j*2*d/lambdas)
    denom=1 + R01(lambdas)*R12(lambdas)*np.exp(2*1.48*np.pi*j*2*d/lambdas)
    return 1+num/denom

#Definition utilisée précédemment
#def RA(lambdas,d):
#    #In the article R coefficient are taken positive. Here they are negative
#    j=complex(0,1)
#    num=(1-R01(lambdas))**2*np.exp(2*1.48*np.pi*j*2*d/lambdas)
#    denom=1 + R01(lambdas)*R12(lambdas)*np.exp(2*1.48*np.pi*j*2*d/lambdas)
#    return 1-num/denom

##Definition prenant en compte la 1ere réflexion entre le tube et la surface d'oxyde
#def RA(lambdas,d):
#    #In the article R coefficient are taken positive. Here they are negative
#    j=complex(0,1)
#    num=R01(lambdas)+R12(lambdas)*np.exp(2*1.48*np.pi*j*2*d/lambdas)
#    denom=1 + R01(lambdas)*R12(lambdas)*np.exp(2*1.48*np.pi*j*2*d/lambdas)
#    return 1+num/denom

def Kramers(En,Imag,step):
    KK=[]
    y=interp1d(En,Imag,kind='cubic') #interp1d is deprecated. The function UnivariateSpline seems to be prefered
                                    #But it is UGLY !! k=5 means cubic interpolation
    #Range only takes integer so we trick it by divinding the max by the value of the step
    #and then multiplying the value of x (or xx) by the step, back
    n=0
    for x in range(int(1.5/step),int(3/step)): # Only compute the values we will plot to save some time
        K=0
        #We don't take the case where x=xx to avoid division by 0 
        for xx in range(1,int(5/step)): #We don't start at 0 to avoid some complications
            if x != xx:
                if xx*step > 0.55 and xx*step < 0.65:
                    A=5
                elif xx*step > 1.05 and xx*step < 1.15:
                    A=3
                elif xx*step > 1.5 and xx*step < 3:
                    A=float(y(xx*step))
                elif xx*step>3.1 and xx*step<3.2:
                    A=0.1
                else :
                    A=0.05
                K+=((xx*step)*A*step)/((xx*step)**2-(x*step)**2)
        K=K*2/np.pi
        KK.append(K)
        n+=1
        printProgressBar(n + 1, 1.5/step, prefix = 'Kramers calculation in Progress:', suffix = 'Complete', length = 50)
    return KK

def inverse(x,x0,A):
    return A/(x-x0)


class Function():
    """Generic class to define the functions
    Other function will inherite from this class"""
    
    def __init__(self,**kwargs):
        self.Initarg = {}
        for k,v in kwargs.items():
            self.Initarg[k] = v
        self.reinit()
        #self.param = [*args]
        self.reprname = "Generic function Class "
        self.name = "No Function : "
    
    def function(self, *args):
        print("this is a generic function")
        print("you need to define a nex function")
        return
    
    def get_values(self):
        return self.param
    
    def reinit(self):
        '''Reinitialize the values of the parameter to the original called ones
        returns self'''
        self.param=[i for i in self.Initarg.values()]
        self.paramnames = [i for i in self.Initarg.keys()]
        return self
    
    def draw(self,x,*args,**kwargs):
        return plt.plot(x,self.function(x,*self.param),*args,**kwargs)
    
    def __repr__(self):
        string=self.reprname + " Params={}".format(str(self.param))
        return string
    
    def __str__(self):
        string = self.name + " Parameters are  " + ", ".join(["=".join([i,str(j)]) for i,j in zip(self.paramnames,self.param)])
        return string
    
    @classmethod
    def new(cls,*args):
        return cls(*args)


class Contrast(Function):
    
    def __init__(self,delta0,ReChi,ImChi,f):
        super().__init__(delta0=delta0,ReChi=ReChi,ImChi=ImChi,f=f)
        #fy est fixé en paramètre dans le cas du modèle complet, sinon on ne converge pas
        self.j=complex(0,1)
        self.lambdas=560
        self.d=100
        self.RE=RE(self.lambdas,self.d)
        self.RA=RA(self.lambdas,self.d)
        self.fy=0.001
        self.name = 'Contrast function'
        self.reprname = self.name
        
    def function(self,delta,delta0,ReChi,ImChi,f):
        r""" Fonction Contrast for the article
        
        .. math:: \frac{|R_A|^4|A|^2}{4|R_E|^2(|\delta-\delta_0+if|^2+f_y^2)}+ \frac{Re(e^{i\phi}A^*R_A^2R_E^2(\delta-\delta_0+if)}{|R_E|^2(|\delta-\delta_0+if|^2+f_y^2))}
        
        """
        phi=1/2
        Chi = ReChi + self.j*ImChi
        num1=np.abs(self.RA)**4*np.abs(Chi)**2
        denom1=4*np.abs(self.RE)**2*(np.abs(delta-delta0+self.j*f)**2+self.fy**2)
        num2=np.real(np.exp(self.j*phi*np.pi)*np.conj(Chi*self.RA**2)*self.RE*(delta-delta0+self.j*f))
        denom2=np.abs(self.RE)**2*(np.abs(delta-delta0+self.j*f)**2+self.fy**2)
        return num1/denom1+num2/denom2

    def function_scat(self,delta,delta0,ReChi,ImChi,f):
        r""" Fonction Contrast for the article
        
        .. math:: \frac{|R_A|^4|A|^2}{4|R_E|^2(|\delta-\delta_0+if|^2+f_y^2)}+ \frac{Re(e^{i\phi}A^*R_A^2R_E^2(\delta-\delta_0+if)}{|R_E|^2(|\delta-\delta_0+if|^2+f_y^2))}
        
        """
        Chi = ReChi + self.j*ImChi
        num1=np.abs(self.RA)**4*np.abs(Chi)**2
        denom1=4*np.abs(self.RE)**2*(np.abs(delta-delta0+self.j*f)**2+self.fy**2)
        return num1/denom1
    
    def draw_scat(self,x,*args,**kwargs):
        return plt.plot(x,self.function_scat(x,*self.param),*args,**kwargs)
    
class Contrast_ideal(Function):    
    def __init__(self,ReChi,ImChi):
        super().__init__(ReChi=ReChi,ImChi=ImChi)
        self.j=complex(0,1)
        self.lambdas=560
        self.d=100
        self.RE=RE(self.lambdas,self.d)
        self.RA=RA(self.lambdas,self.d)
        self.name = 'Contrast ideal function'
        self.reprname = self.name
        
    def function(self,delta,ReChi,ImChi):
        r""" Fonction Contrast for the article
        
        .. math:: \frac{|R_A|^4|A|^2}{4|R_E|^2(|\delta|^2)}+ \frac{Re(e^{i\phi}A^*R_A^2R_E^2(\delta)}{|R_E|^2(|\delta|^2))}
        
        """
        phi=1/2
        Chi = ReChi + self.j*ImChi
        num1=np.abs(self.RA)**4*np.abs(Chi)**2
        denom1=4*np.abs(self.RE)**2*(np.abs(delta)**2)
        num2=np.real(np.exp(self.j*phi*np.pi)*np.conj(Chi*self.RA**2)*self.RE*(delta))
        denom2=np.abs(self.RE)**2*(np.abs(delta)**2)
        return num1/denom1+num2/denom2

class Contrast_NCoherent(Function):
    
    def __init__(self,ReChi,ImChi,fy):
        super().__init__(ReChi=ReChi,ImChi=ImChi.imag,fy=fy)
        self.j=complex(0,1)
        self.lambdas=560
        self.d=100
        self.RE=RE(self.lambdas,self.d)
        self.RA=RA(self.lambdas,self.d)
        self.name = 'Contrast non coherent function'
        self.reprname = self.name
        
    def function(self,delta,ReChi,ImChi,fy):
        r""" Fonction Contrast for the article
        
        .. math:: \frac{|R_A|^4|A|^2}{4|R_E|^2(|\delta|^2+f_y^2)}+ \frac{Re(e^{i\phi}A^*R_A^2R_E^2(\delta)}{|R_E|^2(|\delta|^2+f_y^2))}
        
        """
        Chi = ReChi + self.j*ImChi
        phi=1/2
        num1=np.abs(self.RA)**4*np.abs(Chi)**2
        denom1=4*np.abs(self.RE)**2*(np.abs(delta)**2+fy**2)
        num2=np.real(np.exp(self.j*phi*np.pi)*np.conj(Chi*self.RA**2)*self.RE*(delta))
        denom2=np.abs(self.RE)**2*(np.abs(delta)**2+fy**2)
        return num1/denom1+num2/denom2

class imagLorentz(Function):
    def __init__(self,lambda0,Gamma,A):
        super().__init__(lambda0=lambda0,Gamma=Gamma,A=A)
        self.name='complex lorentz'
        self.reprname = self.name
    
    def function(self,x,x0,gamma,A):
        j=np.complex(0,1)
        return A/((1/x0-1/x)-j*gamma/2)
    
    
class Sin(Function):
    r""" Sin class to generate a sin function.
    
    """
    def __init__(self,x0,Omega,A):
        super().__init__(x0=x0,Omega=Omega,A=A)
        self.name = 'Sin function'
        self.reprname = self.name
        
    def function(self,x,x0,Omega,A):
        r""" Fonction sinus
        .. math:: A*\sin{\omega*(x-x_0)}
        """
        return A*np.sin(Omega*(x-x0))


class Fano(Function):
    r""" Fano class to generate a fano profil function.
    
    """
    def __init__(self,x0,Gamma,A,q=10):
        super().__init__(x0=x0,Gamma=Gamma,A=A,q=q)
        self.name = 'Fano function'
        self.reprname = self.name
        
    def function(self,x,x0,Gamma,A,q):
        r""" Fonction Fano
        .. math:: \frac{A(q.\Gamma+2*(x-x0))}{(\Gamma^2+4*(x-x0)^2)}
        """
        num=q*Gamma + 2*(x-x0)
        denom=Gamma**2 + 4*(x-x0)**2
        return A/q**2*(num**2/denom)
    
    
class Gaussian(Function):
    r""" non-normalized Gaussian class to generate a gaussian function.
    
    """
    def __init__(self,x0,Gamma,A):
        super().__init__(x0=x0,Gamma=Gamma,A=A)
        self.name = 'Gaussian function'
        self.reprname = self.name
        
    def function(self,x,x0,Gamma,A):
        r""" Fonction Gaussian
        
        .. math:: A*e^{\frac{-(x-x_0)^2}{2*\Gamma^2}}
        
        """
        expo=-(x-x0)**2/(2*Gamma**2)
        norm=A
        return norm*np.exp(expo)
    

class GaussianNorm():
    r""" non-normalized Gaussian class to generate a gaussian function.
    
    """
    def __init__(self,x0,Gamma,A):
        super().__init__(x0=x0,Gamma=Gamma,A=A)
        self.name = 'Norm. Gaussian function'
        self.reprname = self.name
        
    def function(self,x,x0,Gamma,A):
        r""" Fonction Gaussian
        
        .. math:: \frac{A}{\Gamma*\sqrt{2}*\pi}*e^{\frac{-(x-x_0)^2}{2*\Gamma^2}}
        
        """
        expo=-(x-x0)**2/(2*Gamma**2)
        norm=A/(Gamma*np.sqrt(2)*np.pi)
        return norm*np.exp(expo)

class Lorentz(Function):
    r""" Lorentz class to generate a lorentzian function.
    
    """
    def __init__(self,x0,Gamma,A):
        super().__init__(x0=x0,Gamma=Gamma,A=A)
        self.name = 'Lorentz function'
        self.reprname = self.name
        
    def function(self,x,x0,Gamma,A):
        r""" Fonction Lorentz
        
        .. math:: a+bx+cx^2+ \frac{A}{1+4*(x-x0)^2/\Gamma^2}
        
        """
        num=A
        denom=1+4*((x-x0)/Gamma)**2
        return num/denom

class Polynom2(Function):
    r""" Polynom of degree 2"""
    def __init__(self,A,x0,y0):
        super().__init__(A=A,x0=x0,y0=y0)
        self.name='Polynom of degree 2'
        self.repr=self.name
        
    def function(self,x,A,x0,y0):
        return A*((x-x0)**2+y0)

class Polynom(Function):
    r""" Polynom class to generate a polynomial function.
    
    """
    def __init__(self,n=1,array=[1,1]):
        self.n=n
        self.Initarg = {'n':self.n}
        self.array=array
        for i,x in enumerate(self.array):
            self.Initarg['x'+str(self.n-i)] = x
        self.reinit()
        super().__init__(**self.Initarg)
        self.name="Polynom function"
        self.reprname=self.name
    
    
    
    def function(self,x,*args):
        r""" Fonction polynom of degree n
        """
        return np.poly1d(args[1:])(x)
    
    def reinit(self):
        '''Reinitialize the values of the parameter to the original called ones
        returns self'''
        self.Initarg = {'n':self.n}
        if len(self.array)<self.n+1:
            while len(self.array)!=self.n+1: self.array.append(1) 
        elif len(self.array)>self.n+1:
            self.array=self.array[:self.n+1]
        for i,x in enumerate(self.array):
            self.Initarg['x'+str(self.n-i)] = x
        self.param=[i for i in self.Initarg.values()]
        self.paramnames = [i for i in self.Initarg.keys()]
        return self
    
    def __str__(self):
        pol=[]
        for i in range(self.n+1):
            pol.append('{}*x^{}'.format(self.param[i+1],self.n-i))
        string=self.name + " of degree {}. Polynom ={}".format(self.n,' + '.join(pol))
        return string


        

def functionsomme(x,*args,**kwargs):
    start=0
    end=0
    fitmul=0
    for func in kwargs['function']:
        end+=len(func.param)
        fitmul = fitmul + func.function(x,*args[start:end])
        start=end
    return fitmul 
        
def fit(X_array,Y_array,bornes=None,function=[Fano(1,1,1,1)]):
    bestv=[]
    fit=[]
    if bornes is None:
        bornes=np.arange(len(X_array))
    Emin=min(X_array[bornes[0]],X_array[bornes[-1]])
    Emax=max(X_array[bornes[0]],X_array[bornes[-1]])
    step=(Emax-Emin)/500
    seed=[]
    funcfit=[]
    for func in function:
        for param in func.param:
            seed.append(param)
        funcfit.append(func)
    #print(seed)
    def fitfit(x,*args):
        return functionsomme(x,*args,function=funcfit)
    a,b=curve_fit(fitfit,X_array[bornes],Y_array[bornes],seed,maxfev=500000)
    bestv=[a]
    start=0
    end=0
    for func in function:
        end+=len(func.param)
        func.param=bestv[0][start:end]
        start=end
    fit.append(fitfit(np.arange(Emin,Emax,step),*a))
    bestv=np.asarray(bestv).reshape(1,len(seed))
    return bestv,fit
    

def fit_fig6(wavelength,data,Angles,function):
    Fit=[[function]]
    n=0
    for lambdas,contrast in zip(wavelength,data.T):
        d=100
        function=[Fit[-1][0].new(*Fit[-1][0].param)]
        function[0].d=d
        function[0].lambdas=lambdas
        function[0].RE=RE(function[0].lambdas,function[0].d)
        function[0].RA=RA(function[0].lambdas,function[0].d)
        a = fit(Angles,contrast,function=function)
        Fit.append(function)
        printProgressBar(n + 1, len(wavelength), prefix = 'Fit in Progress:', suffix = 'Complete', length = 50)
        n+=1
    return Fit

def plot_fit(FREQ,function,sig,table=None,bornes=None,delta=0.003,shift=True,initfunc=True):
    r""" Fit un ensemble de courbe. (optionel trace les courbes fittées dans une figure)
    
    Parameters
    ------------
    FREQ : array_like
        abscisse pour le fit (e.g fit en energie à un angle donnée FREQ = energie; fit en angle a une energie donnée FREQ = angles)
    function : array of object from class functions
        ex : [Fano(800,10,10,10),Lorentz(10,20,30)]
    sig : array
        données à fitter. 
        Typiquement le diagramme de bande complet. La selection des courbes à fitter se fait dans la variable "table"
    table : array
        tableau contenant l'ensemble des numéros des courbes à fitter 
        (ex: np.arange(30,50) va fitter les courbes à partir de la 31eme jusqu'à la 51eme de l'ensemble donné dans sig)
    bornes : array
        tableau contenant les N points sur lequel le fit va être effectué de la forme np.arange(Pixel1,PixelN)
    seed : array
        tableau contenant les seed pour le premier fit. Les suivants sont réalisés de proches en proches en se resservant de la sortie du fit précédant comme seed.
    delta : float, optional
        indique l'écart entre la valeur de x0 trouvé pour le fit et les bornes pour el fit suivant ([x0-delta;x0+delta])
    shift : {True,False}, optional
        True : shift les bornes entre chaque fit; False : garde les bornes initiales entre chaque fit
    
    
    
    Returns
    ---------
    table : array
        identique à l'entrée
    bestval : array
        tableau contenant l'ensemble des paramètres de fit obtenus pour la série de courbe
    Q : array
        tableau contenant l'ensemble des facteurs de qualité obtenus pour la série de courbe
    yfit : array
        tableau contenant l'ensemble des courbes de fit
    """
    bestval=[] 
    yfit=[]
    Q=[]
    if table is None :
        table = np.arange(len(sig))
    if bornes is None:
        bornes=np.arange(len(FREQ))
    print(len(sig),'bornes= ',FREQ[bornes[0]],' ',FREQ[bornes[-1]])
    if initfunc:
        for func in function:
            func.reinit()
    for n,i in enumerate(table):
        a,b=fit(FREQ,sig[i],bornes,function=function)
        bestval.append(a)
        yfit.append(b)
        Q.append(abs(bestval[n][0][0]/bestval[n][0][1]))
        Emin=min(FREQ[bornes[0]],FREQ[bornes[-1]])
        Emax=max(FREQ[bornes[0]],FREQ[bornes[-1]])
        #print(seed[0])
        if shift == True:
            Em=(np.abs(FREQ-a[0][0]-delta)).argmin()
            #print(Em)
            EM=(np.abs(FREQ-a[0][0]+delta)).argmin()
            #print(EM)
            bornes=np.arange(min(Em,EM),max(Em,EM))
        printProgressBar(n + 1, len(table), prefix = 'Fit in Progress:', suffix = 'Complete', length = 50)
    bestval=np.asarray(bestval).reshape(-1,len(a[0]))
    return table,bestval,Q,yfit

def main():
    return

if __name__ == '__main__':
    main()