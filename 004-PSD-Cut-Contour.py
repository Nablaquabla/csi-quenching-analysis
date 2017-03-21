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
textFS = 10
bgColor = 'w'
eColor='k'
colors = ez.colors()


#=============================================================================#
#                                                                             #
#                           Main program starts here                          #
#                                                                             #
#=============================================================================#
mDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Final-Data-Set/Experimental/PSD-Calibration/'
plotDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Plots/'

def upper_cut_func(x):
    return 0.65*(1.0 - np.exp(-(x+0)/0.1))

def lower_cut_func(x):
    return 0.55*(1.0 - np.exp(-(x-0.05)/0.1))

xp = np.linspace(0,1,1000)
# =============================================================================
#    Cf and Cs in one plot
# =============================================================================

acceptance = {}
plt.figure(figsize=(6.5,3.25),edgecolor='k',facecolor='w')

gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
gs.update(left=0.11, right=0.97, hspace=0.00,wspace=0.02,bottom=0.135,top=0.97)

ax = []
ax.append(plt.subplot(gs[0]))
ax.append(plt.subplot(gs[1]))

plt.sca(ax[0])
dataset = 'Cf-Plastic/'
onset,full,tail = np.loadtxt(mDir + dataset + 'EJ-0000',unpack=True)

en = full
psd = (full-tail)/full
nF,bF = np.histogram(en,100,[0,1])
cut = (psd <= upper_cut_func(en)) * (psd >= lower_cut_func(en))
nC,bC = np.histogram(en[cut],100,[0,1])
acceptance['Cf'] = np.array([nC,nF])

plt.scatter(en,psd,marker='+',color=colors.red,alpha=0.1)
l1 = plt.scatter([0],[0],marker='+',color=colors.red,alpha=1,label='$^{252}$Cf')

dataset = 'Cs-Plastic/'
onset,full,tail = np.loadtxt(mDir + dataset + 'EJ-0000',unpack=True)

en = full
psd = (full-tail)/full
nF,bF = np.histogram(en,100,[0,1])
cut = (psd <= upper_cut_func(en)) * (psd >= lower_cut_func(en))
nC,bC = np.histogram(en[cut],100,[0,1])
acceptance['Cs'] = np.array([nC,nF])

plt.scatter(en,psd,marker='+',color=colors.black,alpha=0.1)
l2 = plt.scatter([0],[0],marker='+',color=colors.black,alpha=1,label='$^{137}$Cs')

leg1 = plt.legend(handles=[l1,l2],loc='upper right',framealpha=0,fontsize=legendFS,ncol=2)
vp = leg1._legend_box._children[-1]._children[0]
for c in vp._children:
    c._children.reverse()
vp.align="right"

axLeg = plt.gca().add_artist(leg1)


plotCol = '#5088B2'
l3, = plt.plot(xp,upper_cut_func(xp),c=plotCol,lw=1.5,label='PSD slection window')
l3, = plt.plot(xp,lower_cut_func(xp),c=plotCol,lw=1.5,label='PSD slection window')
plt.fill_between(xp,lower_cut_func(xp),upper_cut_func(xp),color=plotCol,alpha=0.25)

leg2 = plt.legend(handles=[l3],loc='lower right',framealpha=0,fontsize=legendFS)
vp = leg1._legend_box._children[-1]._children[0]
for c in vp._children:
    c._children.reverse()
vp.align="right"

plt.xlim(0,1)
plt.ylim(0.35,1)

plt.ylabel(r'PSD metric $f_\mathrm{psd}$',labelpad=8)
plt.tick_params(axis='x',labelbottom='off')
plt.sca(ax[1])
x = np.linspace(0,1,100)
y = 1.0*acceptance['Cf'][0]/acceptance['Cf'][1]
yerr = np.sqrt(1.0/np.sqrt(acceptance['Cf'][0]) + 1.0/np.sqrt(acceptance['Cf'][1]))*y
plt.errorbar(x,y,yerr=yerr,linestyle='None',marker='.',c=colors.red)

y = 1.0*acceptance['Cs'][0]/acceptance['Cs'][1]
yerr = np.sqrt(1.0/np.sqrt(acceptance['Cs'][0]) + 1.0/np.sqrt(acceptance['Cs'][1]))*y
plt.errorbar(x,y,yerr=yerr,linestyle='None',marker='.',c=colors.black)

plt.axhline(1.0,color='k',ls='dashed')
plt.axhline(0.01,color='k',ls='dashed')
plt.yscale('log',nonposy='clip')
plt.ylim(0.0001,5)
plt.xlabel(r'Current [a.u.]')
plt.ylabel('Acceptance')
plt.savefig(plotDir + 'Plastic-PSD-Cs-And-Cf',dpi=500)
plt.show()





















