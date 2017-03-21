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
import h5py
from scipy.stats import poisson
import colormaps as cmaps
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatter

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
legendFS = 10
labelFS = 10

bgColor = 'w'
eColor='k'
colors = ez.colors()
colorMap = cmaps.viridis

#=============================================================================#
#                                                                             #
#                           Main program starts here                          #
#                                                                             #
#=============================================================================#
dataOutputDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Final-Data-Set/Results/'
plotDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Plots/'

angles = np.asarray([18,21,24,27,33,39,45])
nEnergies = np.asarray([2.87,4.04,4.80,6.18,9.07,12.61,16.37])
nEnergiesSigma = np.asarray([0.56,0.72,0.80,0.90,1.21,1.43,1.81])# * 2.355

qf = []
qfErrPosFit = []
qfErrNegFit = []
qfErrPosEJ = []
qfErrNegEJ = []
qfErrNegTP = []
qfErrPosTP = []

fitType = 'unrestricted'
#fitType = 'restricted'

for a in angles:
    _d = np.loadtxt(dataOutputDir + '%d-%s.dat'%(a,fitType),skiprows=1).T
    qf.append(_d[1][2])
    qfErrNegFit.append(_d[1][2] - _d[2][2])
    qfErrPosFit.append(_d[3][2] - _d[1][2])
    qfErrNegEJ.append(_d[1][2] - np.min(_d[1]))
    qfErrPosEJ.append(np.max(_d[1]) - _d[1][2])
    qfErrNegTP.append(_d[1][2] - np.min(_d[1::6]))
    qfErrPosTP.append(np.max(_d[1::6]) - _d[1][2])

qfErrPosFit = np.asarray(qfErrPosFit)
qfErrNegFit = np.asarray(qfErrNegFit)
qfErrPosEJ = np.asarray(qfErrPosEJ)
qfErrNegEJ = np.asarray(qfErrNegEJ)
qfErrPosTP = np.asarray(qfErrPosTP)
qfErrNegTP = np.asarray(qfErrNegTP)

qf = np.asarray(qf)
qfErrPosTotal = np.sqrt(qfErrPosFit*qfErrPosFit + qfErrPosEJ*qfErrPosEJ + qfErrPosTP*qfErrPosTP + 0.025*0.025*qf*qf)*100.0
qfErrNegTotal = np.sqrt(qfErrNegFit*qfErrNegFit + qfErrNegEJ*qfErrNegEJ + qfErrNegTP*qfErrNegTP + 0.025*0.025*qf*qf)*100.0
qf = qf*100.0

plt.figure(figsize=(6.5,4),edgecolor='k',facecolor='w')
plt.errorbar(nEnergies,qf,yerr=[qfErrNegTotal,qfErrPosTotal],xerr=nEnergiesSigma,marker='.',capsize=0,color=colors.blue,linestyle='None')

#Park 1
#d = np.loadtxt('./QF-Park-CsI-01.txt',skiprows=1).T
#cCol = colors.brown
#plt.errorbar(d[0],d[1],xerr=d[2],yerr=d[3]-d[1],ecolor=cCol,color=cCol,linestyle='None',marker='o',label='Park $et\,al.$, arXiv:nucl-ex/0202014',capsize=0)

#Park 2
#d = np.loadtxt('./QF-Park-CsI-02.txt',skiprows=1).T
#cCol = colors.brown
#plt.errorbar(d[0],d[1],xerr=d[2],yerr=d[3]-d[1],ecolor=cCol,color=cCol,linestyle='None',marker='o')

#Guo
#d = np.loadtxt('./QF-Guo-CsI.txt',skiprows=1).T
#cCol = colors.green
#plt.errorbar(d[0],d[1],xerr=d[2],yerr=d[3]-d[1],ecolor=cCol,color=cCol,linestyle='None',marker='o',label='Guo $et\,al.$, arXiv:1602.04923',capsize=0)

#Juan
#d = np.loadtxt('./QF-Juan-CsI.txt',skiprows=1).T
#cCol = colors.red
#plt.errorbar(d[0],d[1],xerr=d[2],yerr=d[3]-d[1],ecolor=cCol,color=cCol,linestyle='None',marker='o',label='Collar $et\,al.$, arXiv:1407.7524',capsize=0)


#plt.errorbar(nEnergies,qf,yerr=[qfErrNegFit,qfErrPosFit],xerr=nEnergiesSigma,marker='o',capsize=0,color='k',linestyle='None')
#plt.errorbar(nEnergies+0.25,qf,yerr=[qfErrNegEJ,qfErrPosEJ],xerr=nEnergiesSigma,marker='o',capsize=0,color='r',linestyle='None')
#plt.errorbar(nEnergies+0.5,qf,yerr=[qfErrNegTP,qfErrPosTP],xerr=nEnergiesSigma,marker='o',capsize=0,color='g',linestyle='None')
#plt.ylabel(r'Quenching factor w.r.t. 57.6 keV$_\mathrm{ee}$')

plt.legend(loc='upper right',framealpha=0,fontsize=legendFS)
plt.ylabel(r'Scintillation efficiency w.r.t. 57.6 keV$_\mathrm{ee}$ [%]' )
plt.xlabel(r'Recoil energy [keV$_\mathrm{nr}$]')
#plt.xlim(0,90)
#plt.ylim(0,20)
plt.xlim(0,20)
plt.ylim(2,8)
plt.tight_layout(pad=0.25)



flux = []
fluxErrPosFit = []
fluxErrNegFit = []
fluxErrPosEJ = []
fluxErrNegEJ = []
fluxErrNegTP = []
fluxErrPosTP = []
for a in angles:
    _d = np.loadtxt(dataOutputDir + '%d-%s.dat'%(a,fitType),skiprows=1).T
    flux.append(_d[4][2])
    fluxErrNegFit.append(_d[4][2] - _d[5][2])
    fluxErrPosFit.append(_d[6][2] - _d[4][2])
    fluxErrNegEJ.append(_d[4][2] - np.min(_d[4]))
    fluxErrPosEJ.append(np.max(_d[4]) - _d[4][2])
    fluxErrNegTP.append(_d[4][2] - np.min(_d[4::6]))
    fluxErrPosTP.append(np.max(_d[4::6]) - _d[4][2])
fluxErrPosFit = np.asarray(fluxErrPosFit)
fluxErrNegFit = np.asarray(fluxErrNegFit)
fluxErrPosEJ = np.asarray(fluxErrPosEJ)
fluxErrNegEJ = np.asarray(fluxErrNegEJ)
fluxErrPosTP = np.asarray(fluxErrPosTP)
fluxErrNegTP = np.asarray(fluxErrNegTP)
flux = np.asarray(flux)
fluxErrPosTotal = np.sqrt(fluxErrPosFit*fluxErrPosFit + fluxErrPosEJ*fluxErrPosEJ + fluxErrPosTP*fluxErrPosTP)
fluxErrNegTotal = np.sqrt(fluxErrNegFit*fluxErrNegFit + fluxErrNegEJ*fluxErrNegEJ + fluxErrNegTP*fluxErrNegTP)

print flux[-1], fluxErrPosTotal[-1],fluxErrNegTotal[-1]
plt.figure(figsize=(6.5,4),edgecolor='k',facecolor='w')
plt.errorbar(angles,flux/1000.0,yerr=[fluxErrNegFit/1000.0,fluxErrPosFit/1000.0],marker='o',capsize=0,color=colors.blue,linestyle='None')
plt.ylabel('Neutron flux [1000 n /s]')
plt.xlabel('Angle [degree]')
plt.xlim(15,50)
plt.tight_layout(pad=0.25)
plt.show()











