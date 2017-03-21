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
colorMap = cmaps.viridis

#=============================================================================#
#               Poisson smears data using precalculated kernels
#=============================================================================#
# Prepare kernels for poisson smearing
xKernel = np.arange(41)
kernels = np.array([poisson.pmf(xKernel, xx) for xx in xKernel])

# Smear function
def poissonSmear(hist,kernels):
    _t = np.zeros(len(hist))
    for i,h in enumerate(hist):
        _t = _t + kernels[i] * h
    return _t

#=============================================================================#
#                                                                             #
#                           Main program starts here                          #
#                                                                             #
#=============================================================================#
expInputDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Final-Data-Set/Experimental/Residual-Spectra/'
simInputDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Final-Data-Set/Simulations/Output/'
plotDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Plots/'
plotFigures = True

# Choose colors used for plotting both types
expColor = colors.red
simColor = colors.black

# Number of bins used during fitting
nQF = 100
nS = 100

# Calibration data from Am calibration and in-situ SPEQ measurement
speq = {'gauss': {'best': 0.0089, 'negative': 0.0086, 'positive': 0.0092},
        'polya': {'best': 0.0118, 'negative': 0.0115, 'positive': 0.0121}}

runtimes = {45: 165,
            39: 182,
            33: 225,
            27: 207,
            24: 178,
            21: 203,
            18: 71}

lightYield = 13.16

angle = 24
ejThreshold = 450
speType = 'polya'
qType = 'best'

# Read simulated and experimental data
fExp = h5py.File(expInputDir + 'Residual-%dDegree-%dkeV-%s-%s.h5'%(angle,ejThreshold,speType,qType),'r')
fSim = h5py.File(simInputDir + '%d-Sim-Results.h5'%angle,'r')
nExp = fExp['n'][...]
nExpErr = fExp['nerr'][...]
nExpErr[nExpErr == 0] = 1
nExp = np.sum(np.reshape(nExp[:-1],(20,2)),axis=1)
nExpErr = np.reshape(nExpErr[:-1],(20,2))
nExpErr = np.sqrt(np.sum(nExpErr * nExpErr,axis=1))

# Get data that passes energy cut on EJ
ejCut = fSim['EJ-Edep'][...] > ejThreshold
csiEdep = np.sum(fSim['CsI-Edep'][...][ejCut,:],axis=1)

# Get scaling and quenching factor arrays that will be scanned by fitting
scalingConversion = 9e8/runtimes[angle]/60.0
sArr = np.linspace(500,2500,nS)
qfSamplingArr = np.linspace(0.04,0.10,nQF)

# Fit simulation to data
# 1. For each QF in the array calculate the number of PE for each simulated event
# that are being detected
# 2. Bin data by NPE
# 3. For each scaling point calculate the chi square
x2Arr = []
qfArr = []
scArr = []
for qf in qfSamplingArr:
    _npe =  csiEdep* lightYield * qf
    nSim,b = np.histogram(_npe,41,[-0.5,40.5])
    nSim = poissonSmear(nSim,kernels)
    nSim = np.sum(np.reshape(nSim[:-1],(20,2)),axis=1)
    for flux in sArr:
        f = flux / scalingConversion
        x2 = np.sum(((nExp-f*nSim)*(nExp-f*nSim)/(nExpErr**2)))
#        x2 = np.sum((nExp-f*nSim)*(nExp-f*nSim))/np.sum(nExpErr)
#        x2 = np.sum((nExp-f*nSim)*(nExp-f*nSim))
        x2Arr.append(x2)
        qfArr.append(qf)
        scArr.append(flux)

x2Arr = np.asarray(x2Arr)
qfArr = np.asarray(qfArr)
scArr =  np.asarray(scArr)

# Determine best fit (QF & scaling) by getting the lowest chi square
bestQF = qfArr[np.argmin(x2Arr)]
bestScaling = scArr[np.argmin(x2Arr)]
x2Min = x2Arr[np.argmin(x2Arr)]
x2Levels = [x2Min + 1, x2Min + 2 ]

# Print n/s and QF from fits
print bestScaling,bestQF

# Plot figures if wished for
if plotFigures:
    xPlot = np.arange(0.5,39,2)
    # 1) 2d distribution of chi square.
    # x = neutron flux // y = quenching factor
#    plt.figure(figsize=(6.5,6.5),edgecolor='k',facecolor='w')
#    extent = [sArr[0],sArr[-1],qfSamplingArr[0],qfSamplingArr[-1]]
#    h = np.reshape(x2Arr,(nQF,nS))
#    plt.imshow(h,aspect='auto',cmap=colorMap,origin='lower',extent=extent)
#    plt.contour(h,x2Levels,extent=extent,linestyles=['solid','dashed'],colors=['w','w'],origin='lower')
#    plt.scatter([bestScaling],[bestQF],marker='.',color='w')
#    plt.xlim(extent[0],extent[1])
#    plt.ylim(extent[2],extent[3])
#    plt.xlabel('Neutron scaling [a.u.]')
#    plt.ylabel('Quenching factor')
#    plt.tight_layout(pad=0.25)

    # 2) Experimental residual spectrum + best fit simulation
    # x = NPE // y = counts
    print len(xPlot), len(nExp), len(nSim)
    print bestScaling, scalingConversion
    plt.figure(figsize=(6.5,6.5),edgecolor='k',facecolor='w')
    _npe =  csiEdep* lightYield * bestQF
    nSim,b = np.histogram(_npe,41,[-0.5,40.5])
    nSim = poissonSmear(nSim,kernels)
    nSim = np.sum(np.reshape(nSim[:-1],(20,2)),axis=1)
    plt.errorbar(xPlot,nExp,yerr=nExpErr,linestyle='None',color=expColor,capsize=0,marker='o')
    plt.plot(xPlot,nSim*bestScaling/scalingConversion,c=simColor)
    plt.axhline(0,color=colors.black,linestyle='dashed')
    plt.ylabel('Counts / photoelectron')
    plt.xlabel('Number of photoelectrons')
    plt.tight_layout(pad=0.25)

    # 3) Combination of (1) & (2)
#    plt.figure(figsize=(6.5,3.25),edgecolor='k',facecolor='w')
#
#    ax0 = plt.subplot(121)
#    extent = [sArr[0]/1000.0,sArr[-1]/1000.0,qfSamplingArr[0],qfSamplingArr[-1]]
#    h = np.reshape(x2Arr,(nQF,nS))
#    plt.imshow(h,aspect='auto',cmap=colorMap,origin='lower',extent=extent)
#    cb = plt.colorbar(orientation='horizontal')
#
#    plt.contour(h,[np.min(h)+1,np.min(h)+2],extent=extent,linestyles=['solid','dashed'],colors=['w','w'],origin='lower')
#    plt.scatter([bestScaling/1000.0],[bestQF],marker='.',color='w')
#    plt.xlim(extent[0],extent[1])
#    plt.ylim(extent[2],extent[3])
#    plt.xlabel('Neutron flux [1000 n/s]')
#    plt.ylabel('Quenching factor')
#
#    ax1 = plt.subplot(122)
#    _npe =  csiEdep* lightYield * bestQF
#    nSim,b = np.histogram(_npe,41,[-0.5,40.5])
#    nSim = poissonSmear(nSim,kernels)
#    nSim = np.sum(np.reshape(nSim[:-1],(20,2)),axis=1)
#    print len(xPlot)
#    print len(nExp)
#    print len(nSim)
#    plt.errorbar(xPlot,nExp,yerr=nExpErr,linestyle='None',color=expColor,capsize=0,marker='o')
#    plt.plot(xPlot,nSim*bestScaling/scalingConversion,c=simColor)
#    plt.axhline(0,color=colors.black,linestyle='dashed')
#    plt.ylabel('Counts / photoelectron')
#    plt.xlabel('Number of photoelectrons')
#    plt.xlim(0,30)
#    plt.tight_layout(pad=0.25)
plt.show()
























