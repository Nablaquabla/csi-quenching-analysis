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
from scipy import integrate, interpolate
from scipy.misc import comb
import matplotlib.gridspec as gridspec
import h5py

#=============================================================================#
#                           Setup  Plot Settings
#=============================================================================#
mpl.rc('font', **{'sans-serif' : 'Times','family' : 'serif'})
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.shadow'] = False
mpl.rcParams['figure.dpi'] = 100
mpl.rcParams['ytick.major.pad'] = 5
mpl.rcParams['xtick.major.pad'] = 5
mpl.rcParams['axes.formatter.limits'] = (-3,3)
plt.rcParams['pdf.fonttype'] = 42
plt.rc('font', **{'size':12})
legendFS = 12
textFS = 12
bgColor = 'w'
eColor='k'
colors = ez.colors()

#=============================================================================#
#                                                                             #
#                           Main program starts here                          #
#                                                                             #
#=============================================================================#
resOutDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Final-Data-Set/Experimental/Residual-Spectra/'
plotDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Plots/'

#def upper_cut_func(x):
#    return 0.65*(1.0 - np.exp(-(x-0.000)/0.1))
#
#def lower_cut_func(x):
#    return 0.55*(1.0 - np.exp(-(x-0.05)/0.1))
def upper_cut_func(x):
    return 0.65*(1.0 - np.exp(-(x-0.000)/0.14))

def lower_cut_func(x):
    return 0.55*(1.0 - np.exp(-(x-0.05)/0.17))

decay_times = {'fast': 527.0, 'slow': 5600.0}
ratio = 0.41
tMax = 3000.0
def cdf(t):
    pf = 1.0/(1+ratio) * (1 - np.exp(-t/decay_times['fast'])) / (1-np.exp(-tMax/decay_times['fast']))
    ps = ratio/(1+ratio) * (1 - np.exp(-t/decay_times['slow'])) / (1-np.exp(-tMax/decay_times['slow']))
    return pf + ps

def arrival(t,n):
    return 1-(1.0 - cdf(t))**n

speTypeArr = ['polya']
chargeArr = ['positive']
energyThArr = [250,300,350,400,450]

signal_onset = -0.0625
#signal_onset = -0.05
bg_offset = 0.025
bg_outside_cut_lo = 0.025
bg_outside_cut_hi = 0.025

plotFigure = False
npe_max = 40

speq = {'gauss': {'best': 0.0089, 'negative': 0.0086, 'positive': 0.0092},
        'polya': {'best': 0.0117, 'negative': 0.0114, 'positive': 0.0120}}

xnpe = np.arange(0,41,1)
#window = np.array([3000, 2874, 1948, 1246, 892, 692, 566, 479, 415, 366, 327, 296, 270, 249, 230, 215, 201, 189, 178, 168, 160, 152, 145, 138, 133, 127, 122, 118, 113, 109, 106, 102, 99, 96, 93, 90, 88, 85, 83, 81, 79])
#window = np.array([1000,576, 264, 172, 127, 101, 84, 71, 62, 55, 50, 45, 41, 38, 35, 33, 31, 29, 27, 26, 24, 23, 22, 21, 20, 19, 19, 18, 17, 17, 16, 16, 15, 15, 14, 14, 13, 13, 13, 12, 12])
#window = np.array([6000,5977, 1152, 653, 462, 359, 294, 249, 216, 191, 171, 155, 141, 130, 120, 112, 105, 99, 93, 88, 83, 79, 76, 72, 69, 66, 64, 61, 59, 57, 55, 53, 52, 50, 48, 47, 46, 44, 43, 42, 41])
window = np.array([30000, 18871, 5977, 2061, 1152, 828, 653, 541, 462, 404, 359, 324, 294, 270, 249, 231, 216, 203, 191, 180, 171, 162, 155, 148, 141, 135, 130, 125, 120, 116, 112, 108, 105, 102, 99, 96, 93, 90, 88, 86, 83])
window = window / 1000.0 + signal_onset
timing_cut_func = interpolate.interp1d(xnpe,window,'linear',fill_value='extrapolate')

#for angle in [33]:
for angle in [18,21,24,27,33,39,45]:
    print 'Analyzing angle %i'%angle
    # ============== Read data file ==============
    outDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Final-Data-Set/Experimental/%s-Degree/'%angle
    data = np.load(outDir + 'Full-Data.npz')

    # ============== Split data ==============
    charge = data['charge_spe']
    npe = data['npe']
    csi_onset = data['csi_onset']
    npe_pt = data['npe_pt']
    overflow = data['overflow']
    plastic_onset = data['plastic_onset']
    long_int = data['long_int']
    short_int = data['short_int']
    csi_onset = (csi_onset-plastic_onset)*2/1000.0

    # ================== Process EJ299 data =============
    long_int[long_int==0] = 1e-9
    psd = (long_int-short_int)/long_int

    # ============== Overflow cut ==============
    no_overflow = (overflow == 0)

    # ============== Prepare pretrace cut ==============
    pt_cut = (npe_pt == 0)

    # ============== Prepare cut on plastic onset ==============
    plastic_cut = plastic_onset > 500

    # =========== Remove overflows, afterglows and wrong onsets ========
    _tempCut = no_overflow * pt_cut * plastic_cut
    charge = charge[_tempCut]
    psd = psd[_tempCut]
    csi_onset = csi_onset[_tempCut]
    long_int = long_int[_tempCut]

    # ========== Linear & Quadratic calibration ==============
    en = long_int*693.95 + 1.13

    # =============== Energy thresholds ==============
    for energyTh in energyThArr:
        en_cut = en > energyTh

        # =============== PSD cut ===============
        psd_cut = (psd <= upper_cut_func(long_int)) * (psd >= lower_cut_func(long_int))

        # =============== SPE charge type ===============
        for speType in speTypeArr:

            # ============= SPE Type Error ===============
            for speCharge in chargeArr:
                q_mean = speq[speType][speCharge]

                # ============== Prepare cut on plastic onset ==============
                npe_cut = (charge/q_mean) > 0.1

                # =========== Determine current dataset ==============
                cut = psd_cut * npe_cut * en_cut
                cut_no_psd = npe_cut

                nBins = npe_max + 1
                rng = [-0.5,npe_max + 0.5]
                n_side = {}
                b_side = {}
                x_side = {}

                onset = csi_onset[cut]
                q = charge[cut]

                wd = window-signal_onset
                wd[wd>(1+signal_onset)]=(1+signal_onset)
                scaling = 1.0/arrival(wd[1:]*1000.0,xnpe[1:])
                scaling = np.concatenate(([1],scaling))
                current_type = 'Signal'
                cut = (onset > signal_onset)  * (onset < timing_cut_func(q/q_mean)) * (onset < 1.0)
                n_side[current_type],b_side[current_type] = np.histogram(q[cut]/q_mean,nBins,rng)

                current_type = 'Low_BG'
                cut = (onset < signal_onset - bg_offset) * (onset > (bg_outside_cut_lo - 1.0))
                n_side[current_type],b_side[current_type] = np.histogram(q[cut]/q_mean,nBins,rng)

                n_cut = n_side['Signal'] - (n_side['Low_BG']) / 0.9125 * wd
                n_cut_err = np.sqrt(n_side['Signal'] + (n_side['Low_BG'])*wd**2/0.9125**2)

                fOut = h5py.File(resOutDir + 'Residual-%dDegree-%dkeV-%s-%s.h5'%(angle,energyTh,speType,speCharge),'w')
                fOut.create_dataset('/n',data=n_cut,dtype=np.float)
                fOut.create_dataset('/nerr',data=n_cut_err,dtype=np.float)
                fOut.close()

                if plotFigure:
                    plt.figure(figsize=(7,7),edgecolor='k',facecolor='w')
                    gs1 = gridspec.GridSpec(3,3)
                    gs1.update(hspace=0.01, wspace=0.01,left=0.11,right=0.97,bottom=0.09,top=0.98)
                    ax = []
                    ax.append(plt.subplot(gs1[1:, :2])) # Main scatter plot
                    ax.append(plt.subplot(gs1[0, :2])) # Top histogram
                    ax.append(plt.subplot(gs1[1:, 2])) # Side histogram
                    ax.append(plt.subplot(gs1[0,2])) #Top Side - PSD

                    bgcolor = colors.black
                    scolor = colors.red

                    # Plot main frame
                    plt.sca(ax[0])
                    plt.fill_between([-1+bg_outside_cut_lo,signal_onset-bg_offset], [0]*2,[npe_max]*2,color=bgcolor,alpha=0.15)
                    plt.scatter(onset,q/q_mean,c=colors.black,marker='+',alpha=0.2)
                    plt.axvline(signal_onset, linestyle='dashed',color=scolor)
                    plt.axvline(-1+bg_outside_cut_lo, linestyle='dashed',color=bgcolor)
                    plt.axvline(signal_onset-bg_offset, linestyle='dashed',color=bgcolor)
                    plt.yticks(np.arange(5,npe_max,5))
                    plt.xticks([-0.75,0,0.75])
                    plt.xlim(-1,1)
                    plt.ylim(0,npe_max)
                    plt.xlabel(r'$\Delta t$ [$\mu$s]')
                    plt.ylabel('Number of photoelectrons')

                    plt.plot(window,xnpe,c=colors.red,ls='dashed')
                    plt.text(0.070,36,r'$\tau_{99}$',color=colors.red)
                    yfill = np.linspace(0,40,250)
                    plt.fill_betweenx(yfill,[signal_onset]*250,timing_cut_func(yfill),color=colors.red,alpha=0.1)

                    plt.sca(ax[3])
                    cTh = (energyTh-1.13)/693.95
                    xPlot = np.linspace(0,1.4,200)
                    plt.scatter(long_int[cut_no_psd],psd[cut_no_psd],marker='+',alpha=0.1,color=colors.black)

                    plotCol=colors.blue
                    plt.plot(xPlot,upper_cut_func(xPlot),c=plotCol,linewidth=1)
                    plt.plot(xPlot,lower_cut_func(xPlot),c=plotCol,label='PSD cut contour',linewidth=1)
                    plt.vlines(cTh,lower_cut_func(cTh),upper_cut_func(cTh),colors='r')
                    plt.fill_between(xPlot,upper_cut_func(xPlot),lower_cut_func(xPlot),color=plotCol,alpha=0.2)
                    plt.legend(loc='upper right',framealpha=0,fontsize=10)
                    plt.xlim(0,1.4)
                    plt.ylim(0.4,1.0)
                    plt.tick_params(labelbottom='off',labelleft='off')

                    # Plot top histogram
                    nBins_top = 100
                    rng_top = [-1,1]
                    plt.sca(ax[1])
                    n_top,b_top = np.histogram(onset,nBins_top,rng_top)
                    x_top = np.diff(b_top)*0.5 + b_top[:-1]
                    plt.step(x_top,n_top,where='mid',color=colors.black)
                    plt.xticks([-0.75,0,0.75])
                    plt.xlim(-1,1)
                    plt.tick_params(labelbottom='off')
                    plt.ylabel(r'Counts / %.2f $\mu$s'%(1.0*np.diff(rng_top)/nBins_top))
                    plt.locator_params(axis='y',nbins=5)

                    # Plot side histogram
                    plt.sca(ax[2])
                    plt.step(n_side['Signal'],np.arange(npe_max+1)-0.5,color=colors.red,label='Signal',where='pre')
                    plt.step((n_side['Low_BG']) / 0.9125 * wd,np.arange(npe_max+1)-0.5,color=colors.black,label='Background',where='pre')
                    plt.errorbar(n_cut,np.arange(npe_max+1),xerr=n_cut_err,color=colors.blue,label='Residual',linestyle='None',marker='.')
                    plt.axvline(0,color=colors.black,linestyle='dashed')
                    plt.yticks(np.arange(5,npe_max,5))
                    plt.xlim(np.min(n_cut-n_cut_err),np.max(n_cut)*2.0)
                    plt.ylim(0,npe_max)
                    plt.locator_params(axis='x',nbins=4)
                    plt.legend(loc='upper right',framealpha=0,fontsize=12)
                    plt.tick_params(labelleft='off')
                    plt.xlabel('Counts / PE')
                    plt.savefig(plotDir + '%i-degree-data-selection.png'%angle,dpi=500)

plt.show()