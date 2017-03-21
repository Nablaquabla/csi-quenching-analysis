# -*- coding: utf-8 -*-
#=============================================================================#
#                           Imports
#=============================================================================#
import numpy as np
import matplotlib as mpl
import matplotlib.pylab as plt
import easyfit as ef
import easyfy as ez

#=============================================================================#
#                           Setup  Plot Settings
#=============================================================================#
mpl.rc('font', **{'sans-serif' : 'Arial','family' : 'serif'})
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.shadow'] = False
mpl.rcParams['figure.dpi'] = 100
mpl.rcParams['ytick.major.pad'] = 5
mpl.rcParams['xtick.major.pad'] = 5
mpl.rcParams['axes.formatter.limits'] = (-3,3)

plt.rcParams['pdf.fonttype'] = 42
plt.rc('font', **{'size':12})
legendFS = 12
fs = 16
bgColor = 'w'
eColor='k'
colors = ez.colors()

# =============================================================================
#       Prepare fit figure
# =============================================================================
a,e,s = np.loadtxt('simulated-mean-recoil-energies.dat',skiprows=1,unpack=True)

plt.figure(figsize=(6.5,3.25),edgecolor='k',facecolor='w')
plt.errorbar(a,e,yerr=s,linestyle='None',marker='o',capsize=0,color=colors.black)
p0 = [1,1,1]
x2,pars,xfit,yfit = ef.fit('poly2',a,e,s,p0)
xPlot = np.linspace(0,50,100)
plt.plot(xPlot,ef.poly(xPlot,pars[0],2),c=colors.blue,ls='dashed')
plt.ylabel(r'Mean CsI Recoil Energy [keV$_\mathrm{nr}$]')
plt.xlabel(r'Angle $\alpha$ [degree]')
plt.ylim(0,20)
plt.xlim(0,50)
plt.tight_layout(pad=0.25)
plt.show()










