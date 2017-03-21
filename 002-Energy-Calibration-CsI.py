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
legendFS = 10
textFS = 9
bgColor = 'w'
eColor='k'
colors = ez.colors()

#=============================================================================#
#                           Fit Functions
#=============================================================================#
def gauss(x,q,s,m):
    return np.exp(-0.5*((x-m*q)/s)**2/m)

def ff(x,p):
    t1 = p[2]*gauss(x,p[0],p[1],1)
    t2 = p[5]*gauss(x,p[3],p[4],1)
    return t1 + t2 + p[6]

#=============================================================================#
#                                                                             #
#                           Main program starts here                          #
#                                                                             #
#=============================================================================#
i = 0
mDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Final-Data-Set/Experimental/Am-Calibration/'
plotDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Plots/'

amq = np.array([])
_t_amq = np.array([])

if 'Am-Charges.npy' not in os.listdir(mDir):
    with open(mDir + 'Am-0000','r') as f:
        for line in f:
            i+=1
            if i%10000==0:
                print 'Dump ', i
                amq = np.concatenate((amq,_t_amq[:]))
                _t_amq = np.array([])
            _t_amq = np.append(_t_amq,float(line.split()[1]))
    np.save(mDir + 'Am-Charges.npy',amq)
else:
    amq = np.load(mDir + 'Am-Charges.npy')

speType = 'Polya'
speCharges = {'Polya': 0.0118, 'Gauss': 0.0089}
speq = speCharges[speType]

print '-------------------------------------'
print '\tAm-241 analysis'
print '-------------------------------------'
nt,bt = np.histogram(amq,750,[0,15])
xt = bt[:-1] + 0.5*np.diff(bt)

plt.figure(figsize=(6.5,3.25),edgecolor='k',facecolor='w')
plt.subplots_adjust(bottom=0.1325, left=0.095, top = 0.970, right=0.99,wspace=0,hspace=0)
plt.subplot(211)
plt.step(xt,nt,where='mid',c=colors.black)
c = (xt>6.75) * (xt<12)
p0 = [9.33,0.3,550,8.0,0.5,80,0]
x2,pars,xfit,yfit = ef.arbFit(ff,xt[c],nt[c],'Poisson',p0)
plt.plot(xfit,yfit,c='r')
plt.plot(xfit,pars[0][2]*gauss(xfit,pars[0][0],pars[0][1],1),c=colors.red,ls='dashed')
print '\nAm-241 peak at: %.3f +/- %.3f a.u.'%(pars[0][0], pars[2][0])
print 'Am-127 peak width: %.3f +/- %.3f a.u.\n'%(pars[0][1], pars[2][1])

plt.plot(xfit,pars[0][5]*gauss(xfit,pars[0][3],pars[0][4],1),c=colors.red,ls='dashed')
print '\nL-shell peak at: %.3f +/- %.3f a.u.'%(pars[0][3], pars[2][3])
print 'L-shell peak width: %.3f +/- %.3f a.u.\n'%(pars[0][4], pars[2][4])

plt.text(4.8,265,'K-shell',va='center',ha='center',color=colors.red,fontsize=textFS)
plt.text(9.4,130,'L-shell',va='center',ha='center',color=colors.red,fontsize=textFS)
plt.text(9.3,660,'59.54 keV',va='center',ha='center',color=colors.red,fontsize=textFS)
c = (xt>3.2) * (xt<6.6)
p0 = [200,4.8,0.2,0]
x2,pars,xfit,yfit = ef.fit('gauss',xt[c],nt[c],'Poisson',p0)
print '\nK-shell peak at: %.3f +/- %.3f a.u.'%(pars[0][1], pars[2][1])
print 'K-shell peak width: %.3f +/- %.3f a.u.\n'%(pars[0][2], pars[2][2])

plt.plot(xfit,yfit,c=colors.red)
plt.ylabel('Counts / 0.02 a.u.')
plt.tick_params(axis='x',labelbottom='off')
plt.yticks([0,200,400,600,800])
plt.xlim(0,15)
plt.ylim(0,800)

print '-------------------------------------'
print '\tI-127 analysis'
print '-------------------------------------'
x = np.array([])
y = np.array([])
angles = [18,21,24,27,33,39,45]
for a in angles:
    outDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Final-Data-Set/Experimental/%s-Degree/'%a
    data = np.load(outDir + 'Full-Data.npz')

    # ============== Read Data ==============
    data = np.load(outDir + 'Full-Data.npz')
    charge = data['charge_spe']
    npe = data['npe']
    csi_onset = data['csi_onset']
    npe_pt = data['npe_pt']
    t50 = data['t50']
    t1090 = data['t1090']
    overflow = data['overflow']
    plastic_onset = data['plastic_onset']
    long_int = data['long_int']
    short_int = data['short_int']
    csi_onset = (csi_onset-plastic_onset)*2 / 1000.0

    # ============== Overflow cut ==============
    no_overflow = (overflow == 0)

    # ============== Prepare pretrace cut ==============
    pt_cut = (npe_pt == 0)

    # ============== Prepare t50 cut ==============
    t50_cut = (t50) > 10

    # ============== Prepare t10_90 cut ==============
    t10_90_cut = (t50) < ((t1090) - 10)

    # ============== Prepare cut on plastic onset ==============
    npe_cut = (np.round(charge/speq) > 0.001)

    # ============== Apply cuts ==============
    b = (t50_cut * t10_90_cut * pt_cut * no_overflow * npe_cut)

    x = np.concatenate((x,csi_onset[b]))
    y = np.concatenate((y,charge[b]))

plt.subplot(212)
nBins = 100
rng = [0,20]
n,b = np.histogram(y,nBins,rng)
x = b[:-1] + 0.5*np.diff(b)
yerr = np.sqrt(n)
yerr[yerr==0] = 1
plt.errorbar(x,n,yerr=yerr,c=colors.black,marker='.',linestyle='None',capsize=0)
c = (x>6.5) * (x<11.5)
p0 = [100,8.9,0.1,0]

x2,pars,xfit,yfit = ef.fit('gauss',x[c],n[c],'Poisson',p0)
print '\nI-127 peak at: %.3f +/- %.3f a.u.'%(pars[0][1], pars[2][1])
print 'I-127 peak width: %.3f +/- %.3f a.u.\n'%(pars[0][2], pars[2][2])

plt.plot(xfit,yfit,c=colors.red)
plt.text(8.85,135,'57.6 keV',va='center',ha='center',color=colors.red,fontsize=textFS)
plt.ylim(0,170)
plt.xlim(0,15.0)
plt.yticks([0,30,60,90,120,150])
plt.xlabel('Current [a.u.]')
plt.ylabel('Counts / 0.2 a.u.')
plt.savefig(plotDir + 'CsI-Energy-Calibration.png',dpi=400)
plt.savefig(plotDir + 'CsI-Energy-Calibration.pdf',dpi=200)
plt.show()





