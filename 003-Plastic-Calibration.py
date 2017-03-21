# -*- coding: utf-8 -*-
#=============================================================================#
#                           Imports
#=============================================================================#
import numpy as np
import matplotlib as mpl
import matplotlib.pylab as plt
import easyfit as ef
import easyfy as ez
import colormaps as cmaps
from scipy import integrate
from scipy.misc import comb
import os
import matplotlib.gridspec as gridspec

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
plt.rc('font', **{'size':10})
legendFS = 9
textFS = 9
bgColor = 'w'
eColor='k'
colors = ez.colors()

#=============================================================================#
#                           Fit Functions
#=============================================================================#
def ff(x,p):
    return p[0]*np.exp(-0.5*((x-p[1])/p[2])**2) + p[3]*x + p[4]


#=============================================================================#
#                                                                             #
#                           Main program starts here                          #
#                                                                             #
#=============================================================================#
charge_scale_factor = 1
mDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Final-Data-Set/Experimental/PSD-Calibration/'
plotDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Plots/'

plt.figure(figsize=(6.5,4.5),edgecolor='k',facecolor='w')
gs = gridspec.GridSpec(2,2)

gs.update(left=0.085, right=0.89, top=0.88, bottom=0.115, wspace=0.02*9/12.0, hspace=0.02)

ax = []

ax.append(plt.subplot(gs[0,0]))
ax.append(plt.subplot(gs[1,0]))
ax.append(plt.subplot(gs[0,1]))
ax.append(plt.subplot(gs[1,1]))

# Plot Cs-137 spectrum
plt.sca(ax[0])
dataset = 'Cs-Plastic/'
onset,full,tail = np.loadtxt(mDir + dataset + 'EJ-0000',unpack=True)
n,b = np.histogram(full,120,[0,1.2])
x = b[:-1] + 0.5*np.diff(b)
plt.step(x*charge_scale_factor,n/1000.0,where='mid',color=colors.black)
plt.text(0.66*charge_scale_factor,3.05,'477 keV',va='center',ha='right',fontsize=textFS,color=colors.black)

c = (x > 0.45) * (x < 0.75)
p0 = [100,0.6,0.05,-100,0]
x2,pars,xfit,yfit = ef.arbFit(ff,x[c],n[c],'Poisson',p0)
print pars[0][1],pars[2][1], pars[0][2]
plt.axvline(pars[0][1]+pars[0][2],c='k',ls='dashed')
plt.plot(xfit*charge_scale_factor ,yfit/1000.0,color=colors.red,lw=2)
plt.yticks([0,1,2,3])
plt.xlim(0,0.9)
plt.ylim(0,3.5)
plt.xticks([0,0.2,0.4,0.6,0.8])
plt.text(0.075,0.275,r'$^{137}$Cs',ha='left',va='center',fontsize=textFS)
plt.tick_params(axis='x',labelbottom='off',labeltop='on')
ax[0].xaxis.set_label_position("top")
plt.ylabel(r'Cts$\times$10$^3$/ 0.01 a.u.')
plt.xlabel('Current [a.u.]',labelpad=10)

# Plot Na-22 spectrum
plt.sca(ax[1])
dataset = 'Na-Plastic/'
onset,full,tail = np.loadtxt(mDir + dataset + 'EJ-0000',unpack=True)
n,b = np.histogram(full,120,[0,1.2])
x = b[:-1] + 0.5*np.diff(b)
plt.step(x*charge_scale_factor,n/1000.0,where='mid',color=colors.black)
plt.text(0.50*charge_scale_factor,3.05,'341 keV',va='center',ha='left',fontsize=textFS,color=colors.black)
plt.text(0.075,0.275,r'$^{22}$Na',ha='left',va='center',fontsize=textFS)
c = (x > 0.29) * (x < 0.6)
p0 = [500,0.41,0.04,-100,0]
x2,pars,xfit,yfit = ef.arbFit(ff,x[c],n[c],'Poisson',p0)
print pars[0][1],pars[2][1], pars[0][2]
plt.axvline(pars[0][1]+pars[0][2],c='k',ls='dashed')
plt.plot(xfit*charge_scale_factor,yfit/1000.0,color=colors.red,lw=2)
plt.ylim(0,3.5)
plt.xlim(0,0.9)
plt.xticks([0,0.2,0.4,0.6,0.8])
plt.yticks([0,1,2,3])
plt.xlabel('Current [a.u.]')
plt.ylabel(r'Cts$\times$10$^3$/ 0.01 a.u.')

plt.sca(ax[2])
xsd = np.loadtxt('hydrogen-xs.dat',skiprows=3).T
xsd[0] = xsd[0]*1000.0
l1, = plt.plot(xsd[0],xsd[2],c=colors.blue,label='Photoelectric')
l2, = plt.plot(xsd[0],xsd[1],c=colors.green,label='Inc. Scattering')
l3, = plt.plot(xsd[0],xsd[3],c=colors.red,label='Pair Prod.')

xsd = np.loadtxt('carbon-xs.dat',skiprows=3).T
xsd[0] = xsd[0]*1000.0
plt.plot(xsd[0],xsd[1],ls='dashed',c=colors.green)
plt.plot(xsd[0],xsd[2],ls='dashed',c=colors.blue)
plt.plot(xsd[0],xsd[3],ls='dashed',c=colors.red)

leg1 = plt.legend(handles=[l1,l2],loc='lower left',framealpha=0,fontsize=legendFS)
axLeg = plt.gca().add_artist(leg1)
leg2 = plt.legend(handles=[l3],loc='lower right', framealpha=0,fontsize=legendFS)

plt.yscale('log',nonposy='mask')
plt.xscale('log',nonposy='mask')
plt.xlabel('Photon energy [keV]')
plt.ylabel('barns / atom',rotation=270,labelpad=11)
plt.yticks([1e-11,1e-8,1e-5,1e-2,1e1])
plt.ylim(1e-13,1e1)
plt.tick_params(axis='both',labelbottom='off',labeltop='on',labelleft='off',labelright='on')
ax[2].yaxis.set_label_position("right")
ax[2].xaxis.set_label_position("top")

# =============================================================================
#     Prepare energy calibration plot
# =============================================================================
en = np.asarray([0.0,341.0,477])
charge = np.asarray([0.0,0.425798+0.0566,0.62208175+0.0723])
charge_err = np.asarray([0.05, 0.056627,0.07262549])

def quadratic_calib(x,p):
    return p[0]*x**2 + p[1]*x

xp = np.linspace(0,700,150)

plt.sca(ax[3])
plt.errorbar(en,charge*charge_scale_factor,yerr=charge_err*charge_scale_factor,marker='o',linestyle='None',color=colors.black,capsize=0)

p0 = [0.01,0]
x2,parsL,xfit,yfit = ef.fit('line',en,charge,charge_err,p0)
plt.plot(xp,ef.line(xp,parsL[0])*charge_scale_factor,color=colors.red,label='linear + offset')
print parsL[0]

p0 = [0.0,0.01]
x2,parsQ,xfit,yfit = ef.arbFit(quadratic_calib,en,charge,charge_err,p0)
plt.plot(xp,quadratic_calib(xp,parsQ[0])*charge_scale_factor,color=colors.green,label='quadratic')
print parsQ[0]


plt.legend(loc='upper left',framealpha=0,fontsize=legendFS)
plt.xlabel('Energy [keV]')
plt.ylabel('Current [a.u.]',rotation=270,labelpad=23)
plt.tick_params(axis='y',labelleft='off',labelright='on')
ax[3].yaxis.set_label_position("right")
plt.xticks([0,200,400,600])
plt.ylim(-0.1,1.1)
plt.savefig(plotDir + 'Plastic-Energy-Calibration-Fit.png',dpi=400)
plt.savefig(plotDir + 'Plastic-Energy-Calibration-Fit.pdf',dpi=200)

plt.show()




























