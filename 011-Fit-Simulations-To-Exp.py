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
#               Poisson smears data using precalculated kernels
#=============================================================================#
# Prepare kernels for poisson smearing
xPlot = np.arange(41)
kernels = np.array([poisson.pmf(xPlot, xx) for xx in xPlot])

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
dataOutputDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Final-Data-Set/Results/'
plotDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Plots/Fits/'
plotFigures = True
saveData = True

# Choose colors used for plotting both types
expColor = colors.red
simColor = colors.black

fileLabel = {True: 'unrestricted',
             False: 'restricted'}

# Number of bins used during fitting
nQF = 400
nS = 400

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

fitRegionArr = {45: [1,30],
                39: [1,30],
                33: [1,25],
                27: [1,20],
                24: [1,18],
                21: [1,15],
                18: [1,12]}

speType = 'polya'
lightYield = 13.16
#freeNeutronFluxScaling = False
freeNeutronFluxScaling = True

qfSamplingArr = np.linspace(0.03,0.10,nQF)

if freeNeutronFluxScaling:
    sArr = np.linspace(500,2500,nS)
else:
#    sArr = np.linspace(1583-np.sqrt(437**2 + 158**2), 1583+np.sqrt(260**2 + 158**2),nS)
    sArr = np.linspace(1583 - 158, 1583 + 158,nS)

angleArr = [18,21,24,27,33,39,45]
#angleArr = [18]
ejThresholdArr = [250,300,350,400,450]
#ejThresholdArr = [450]
qTypeArr = ['best','negative','positive']
#qTypeArr = ['negative']
for angle in angleArr:
    fitData = {'S': {}, 'QF': {}}
    for ejThreshold in ejThresholdArr:
        print angle, ejThreshold
        fitData['S'][ejThreshold] = []
        fitData['QF'][ejThreshold] = []

        for qType in qTypeArr:

            # Read simulated and experimental data
            fExp = h5py.File(expInputDir + 'Residual-%dDegree-%dkeV-%s-%s.h5'%(angle,ejThreshold,speType,qType),'r')
            fSim = h5py.File(simInputDir + '%d-Sim-Results.h5'%angle,'r')
            nExp = fExp['n'][...]
            nExpErr = fExp['nerr'][...]
            nExpErr[nExpErr == 0] = 1

            # Get data that passes energy cut on EJ
            ejCut = fSim['EJ-Edep'][...] > ejThreshold
            csiEdep = np.sum(fSim['CsI-Edep'][...][ejCut,:],axis=1)

            # Get scaling and quenching factor arrays that will be scanned by fitting
            scalingConversion = 9e8/runtimes[angle]/60.0

            # Fit simulation to data
            # 1. For each QF in the array calculate the number of PE for each simulated event
            # that are being detected
            # 2. Bin data by NPE
            # 3. For each scaling point calculate the chi square
            x2Arr = []
            qfArr = []
            scArr = []
            fitRegion = fitRegionArr[angle]
            nExpFit = nExp
            nExpErrFit = nExpErr
#            nExpFit = np.sum(np.reshape(nExp[:-1],(20,2)),axis=1)
#            nExpErrFit = np.reshape(nExpErr[:-1],(20,2))
#            nExpErrFit = np.sqrt(np.sum(nExpErrFit * nExpErrFit,axis=1))
            for qf in qfSamplingArr:
                _npe =  csiEdep* lightYield * qf
                nSim,b = np.histogram(_npe,41,[-0.5,40.5])
                nSimFit = poissonSmear(nSim,kernels)
                nSimFitErr = np.sqrt(nSimFit)
                nSimFitErr[nSimFitErr == 0] = 1
#                nSimFit = np.sum(np.reshape(nSim[:-1],(20,2)),axis=1)
                for flux in sArr:
                    f = flux / scalingConversion
#                    x2 = np.sum(((nExpFit-f*nSimFit)*(nExpFit-f*nSimFit)/(nExpErrFit**2))[fitRegion[0]:fitRegion[1]])
                    x2 = np.sum(((nExpFit-f*nSimFit)*(nExpFit-f*nSimFit)/(nExpErrFit**2 + f*nSimFitErr**2))[fitRegion[0]:fitRegion[1]])
#                    x2 = np.sum(((nExpFit-f*nSimFit)*(nExpFit-f*nSimFit)/(f*nSimFit))[fitRegion[0]:fitRegion[1]])
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
            x2PlotLims = [int(np.ceil(np.min(x2Arr))),int(np.ceil(np.max(x2Arr)))]
            x2Ticks = [x2PlotLims[0],
                      int(np.exp(0.33*(np.log(x2PlotLims[1])-np.log(x2PlotLims[0])) + np.log(x2PlotLims[0]))),
                      int(np.exp(0.66*(np.log(x2PlotLims[1])-np.log(x2PlotLims[0])) + np.log(x2PlotLims[0]))),
                      x2PlotLims[1]]

            # Print n/s and QF from fits
            fitData['S'][ejThreshold].append(bestScaling)
            fitData['QF'][ejThreshold].append(bestQF)

            # Plot figures if wished for
            if plotFigures:
                # Left: 2d distribution of chi square. (x = neutron flux // y = quenching factor)
                # Right: Experimental residual spectrum + best fit simulation (x = NPE // y = counts)
                fig = plt.figure(figsize=(6.5,3.25),edgecolor='k',facecolor='w')
                ax0 = plt.subplot(121)
                extent = [sArr[0]/1000.0,sArr[-1]/1000.0,qfSamplingArr[0],qfSamplingArr[-1]]
                h = np.reshape(x2Arr,(nQF,nS))
                plt.imshow(h,aspect='auto',cmap=colorMap,origin='lower',extent=extent,norm=LogNorm(vmin=x2PlotLims[0],vmax=x2PlotLims[1]))
                plt.text(0.7,0.0945,r'$\chi^2$',ha='center',va='center',fontsize=legendFS)
                formatter = LogFormatter(10, labelOnlyBase=False)
                cbaxes = fig.add_axes([0.17, 0.87, 0.25, 0.050])
                cb=plt.colorbar(cax=cbaxes, orientation='horizontal',ticks=x2Ticks, format=formatter)
                plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), fontsize=labelFS,color=colors.black)

                plt.sca(ax0)
                plt.xlim(extent[0],extent[1])
                plt.ylim(extent[2],extent[3])
                plt.xlabel('Neutron flux [1000 n/s]')
                plt.ylabel('Quenching factor')
                ctr = plt.contour(h,[np.min(h)+1,np.min(h)+2],extent=extent,linestyles=['solid','dotted'],colors=['w','w'],origin='lower')

                # Extract sigma levels from contour line
                p = ctr.collections[0].get_paths()[0]
                v = p.vertices
                x = v[:,0]
                y = v[:,1]
                maxScaling = np.max(x) * 1000.0
                minScaling = np.min(x) * 1000.0
                maxQF = np.max(y)
                minQF = np.min(y)
                fitData['S'][ejThreshold].append(minScaling)
                fitData['S'][ejThreshold].append(maxScaling)
                fitData['QF'][ejThreshold].append(minQF)
                fitData['QF'][ejThreshold].append(maxQF)

                # Plot best fit values
                plt.scatter([bestScaling/1000.0],[bestQF],marker='.',color='w')

                # Plot exp residual and scaled simulated spectrum
                ax1 = plt.subplot(122)
                _npe =  csiEdep* lightYield * bestQF
                nSim,b = np.histogram(_npe,41,[-0.5,40.5])
                nSim = poissonSmear(nSim,kernels)
                plt.errorbar(xPlot,nExp,yerr=nExpErr,linestyle='None',color=expColor,capsize=0,marker='o')
    #            plt.plot(xPlot,nSim*bestScaling/scalingConversion,c=simColor)
                plt.step(xPlot,nSim*bestScaling/scalingConversion,c=simColor,where='mid')
                plt.axhline(0,color=colors.black,linestyle='dashed')
                plt.ylabel('Counts / photoelectron')
                plt.xlabel('Number of photoelectrons')
                plt.xlim(0,30)
                plt.tight_layout(pad=0.25)
                plt.savefig(plotDir + '%d-%d-%s.png'%(angle,ejThreshold,qType),dpi=200)
                plt.savefig(plotDir + '%d-%d-%s.pdf'%(angle,ejThreshold,qType),dpi=200)
                plt.close('all')
    #            plt.show()

    if saveData:
        with open(dataOutputDir + '%s-%s.dat'%(angle,fileLabel[freeNeutronFluxScaling]),'w') as f:
            f.write('EJ-Threshold\tBest-QF\tBest-QF-Err1\tBest-QF-Err2\tBest-S\tBest-S-Err1\tBest-S-Err2\tNeg-QF\tNeg-QF-Err1\tNeg-QF-Err2\tNeg-S\tNeg-S-Err1\tNeg-S-Err2\tPos-QF\tPos-QF-Err1\tPos-QF-Err2\tPos-S\tPos-S-Err1\tPos-S-Err2\n')
            for ejThreshold in ejThresholdArr:
                f.write('%d\t%.5f\t%.5f\t%.5f\t%.1f\t%.1f\t%.1f\t%.5f\t%.5f\t%.5f\t%.1f\t%.1f\t%.1f\t%.5f\t%.5f\t%.5f\t%.1f\t%.1f\t%.1f\n'%(ejThreshold,
                 fitData['QF'][ejThreshold][0],fitData['QF'][ejThreshold][1],fitData['QF'][ejThreshold][2],fitData['S'][ejThreshold][0],fitData['S'][ejThreshold][1],fitData['S'][ejThreshold][2],
                 fitData['QF'][ejThreshold][3],fitData['QF'][ejThreshold][4],fitData['QF'][ejThreshold][5],fitData['S'][ejThreshold][3],fitData['S'][ejThreshold][4],fitData['S'][ejThreshold][5],
                 fitData['QF'][ejThreshold][6],fitData['QF'][ejThreshold][7],fitData['QF'][ejThreshold][8],fitData['S'][ejThreshold][6],fitData['S'][ejThreshold][7],fitData['S'][ejThreshold][8]))





