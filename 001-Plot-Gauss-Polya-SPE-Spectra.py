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
bgColor = 'w'
eColor='k'
colors = ez.colors()

#=============================================================================#
#                           Fit Functions
#=============================================================================#
def polya(x,q,t,m):
    return (1.0*x/q)**(m*(t+1.0)-1.0) * np.exp(-(t+1.0)*x/q) / (np.exp(-(t+1.0)))

def gauss(x,q,s,m):
    return np.exp(-0.5*((x-m*q)/s)**2/m)

def gFit3(x,p):
    _t1 = p[2] * gauss(x,p[0],p[1],1)
    _t2 = p[3] * gauss(x,p[0],p[1],2)
    _t3 = p[4] * gauss(x,p[0],p[1],3)
    _tn = p[5] * np.exp(-x/p[6])
    return (_t1 + _t2 + _t3 + _tn) /( 1.0 + np.exp(- p[7] * (x - p[8]))) + p[9]

def pFit3(x,p):
    _t1 = p[2] * polya(x,p[0],p[1],1)
    _t2 = p[3] * polya(x,p[0],p[1],2)
    _t3 = p[4] * polya(x,p[0],p[1],3)
    _tn = p[5] * np.exp(-x/p[6])
    return (_t1 + _t2 + _t3 + _tn) /( 1.0 + np.exp(- p[7] * (x - p[8]))) + p[9]

#=============================================================================#
#                                                                             #
#                           Main program starts here                          #
#                                                                             #
#=============================================================================#
speq = np.array([])
_t_speq = np.array([])
i = 0
mDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Final-Data-Set/Experimental/Am-Calibration/'
plotDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Plots/'

if 'SPE-Charges.npy' not in os.listdir(mDir):
    with open(mDir + 'SPE-0000','r') as f:
        for line in f:
            i+=1
            if i%10000==0:
                print 'Dump ', i
                speq = np.concatenate((speq,_t_speq[:]))
                _t_speq = np.array([])
            _t_speq = np.concatenate((_t_speq,np.array(line.split(),dtype=float)))
    np.save(mDir + 'SPE-Charges.npy',speq)
else:
    speq = np.load(mDir + 'SPE-Charges.npy')

xmin = 0.0015
xmax = 0.025
downsample = 1
fitPolya = True
fitGauss = True

# =============================================================================
#       Prepare fit figure
# =============================================================================
plt.figure(figsize=(6.5,3.25),edgecolor='k',facecolor='w')

gs = gridspec.GridSpec(3, 2)
gs.update(left=0.1, right=0.99, hspace=0.0,wspace=0,top=0.975,bottom=0.125)

ax = []
ax.append(plt.subplot(gs[:2, 0]))
ax.append(plt.subplot(gs[2, 0]))
ax.append(plt.subplot(gs[:2, 1]))
ax.append(plt.subplot(gs[2, 1]))

plt.sca(ax[0])
n,b = np.histogram(speq,358,[0,0.07])
x = 0.5*np.diff(b) + b[:-1]
n = n[::downsample]
x = x[::downsample]
nerr = np.sqrt(n)
nerr[nerr==0]=1
c = (x>xmin)*(x<xmax)
p0 = [0.0117, 2.6, 3e4, 0, 0, 1e6, 0.5e-3, 5e3, 1e-3, 0]
lims = [[0,1,0,0,0],[1,1,0,0,0],[2,1,0,0,0],[3,1,0,0,0],[4,1,0,0,0],[5,1,0,0,0],[6,1,0,0,0],[7,1,0,0,0],[8,1,0,0,0],[9,1,0,0,0]]
x2,pars,xfit,yfit = ef.arbFit(pFit3,x[c],n[c],nerr[c],p0,lims)
print 'Polya SPE charge: ', pars[0][0]
plt.errorbar(x,n/10000.0,yerr=nerr/10000.0,linestyle='None',color=colors.black,marker='+')
plt.plot(xfit,yfit/10000.0,c=colors.red,label=r'Polya')
plt.plot(xfit,pars[0][5] * np.exp(-xfit/pars[0][6])/10000.0,c=colors.red,ls='dotted')
plt.plot(xfit,pars[0][2]*polya(xfit,pars[0][0],pars[0][1],1)/10000.0,c=colors.red,ls='dashed')
plt.plot(xfit,pars[0][3]*polya(xfit,pars[0][0],pars[0][1],2)/10000.0,c=colors.red,ls='-.')
plt.ylabel(r'Counts$\times10^4$/ 2$\times10^{-4}$ a.u.',labelpad=15)
plt.xlim(0,xmax)
plt.ylim(0,4)
plt.tick_params(axis='x',which='both',labelbottom='off')
plt.legend(loc='upper right',framealpha=0,fontsize=legendFS)
plt.locator_params(axis='y',nticks=3)
plt.yticks([0,1,2,3,4])

plt.sca(ax[1])
plt.errorbar(x[c],(n[c]/pFit3(x[c],pars[0])),yerr=(n[c]/pFit3(x[c],pars[0]))*np.sqrt((nerr[c]/n[c])**2+ (1.0/pFit3(x[c],pars[0]))**2),marker='+',color=colors.red,linestyle='None')
plt.axhline(1,c=colors.black,linestyle='dashed')
plt.xlabel('Current [a.u.]')
plt.ylabel(r'Ratio')
plt.xlim(0,xmax)
plt.ylim(0.85,1.15)
plt.yticks([0.9,1,1.1])
plt.xticks([0,0.010,0.020])

plt.sca(ax[2])
n,b = np.histogram(speq,358,[0,0.07])
x = 0.5*np.diff(b) + b[:-1]
n = n[::downsample]
x = x[::downsample]
nerr = np.sqrt(n)
nerr[nerr==0]=1
c = (x>xmin)*(x<xmax)
p0 = [  9.05748909e-03, 6.00762752e-03, 3.08005539e+04, 1.99900708e+03, 6.06466466e+02, 5.08609824e+07, 2.13064956e-04, 7.70293350e+04, 1.63976066e-03, 2.55175243e+01]
lims = [[0,1,0,0,0],[1,1,0,0,0],[2,1,0,0,0],[3,1,0,0,0],[4,1,0,0,0],[5,1,0,1e4,0],[6,1,0,0,0],[7,1,0,0,0],[8,1,0,0,0],[9,1,0,0,0]]
x2,pars,xfit,yfit = ef.arbFit(gFit3,x[c],n[c],nerr[c],p0,lims)
print 'Gauss SPE charge: ', pars[0][0]
plt.errorbar(x,n/10000.0,yerr=nerr/10000.0,linestyle='None',color=colors.black,marker='+')
plt.plot(xfit,yfit/10000.0,c=colors.green,label=r'Gauss')
plt.plot(xfit,pars[0][5]*np.exp(-xfit/pars[0][6])/10000.0,c=colors.green,ls='dotted')
plt.plot(xfit,pars[0][2]*gauss(xfit,pars[0][0],pars[0][1],1)/10000.0,c=colors.green,ls='dashed')
plt.plot(xfit,pars[0][3]*gauss(xfit,pars[0][0],pars[0][1],2)/10000.0,c=colors.green,ls='-.')
plt.xlim(0,xmax)
plt.ylim(0,4)
plt.tick_params(axis='x',which='both',labelbottom='off')
plt.tick_params(axis='y',which='both',labelleft='off')
plt.legend(loc='upper right',framealpha=0,fontsize=legendFS)
plt.yticks([0,1,2,3,4])

plt.sca(ax[3])
plt.errorbar(x[c],(n[c]/gFit3(x[c],pars[0])),yerr=(n[c]/gFit3(x[c],pars[0]))*np.sqrt((nerr[c]/n[c])**2+ (1.0/gFit3(x[c],pars[0]))**2),marker='+',color=colors.green,linestyle='None')
plt.axhline(1,c=colors.black,linestyle='dashed')
plt.xlabel('Current [a.u.]')
plt.tick_params(axis='y',which='both',labelleft='off')
plt.xlim(0,xmax)
plt.ylim(0.85,1.15)
plt.yticks([0.9,1,1.1])
plt.xticks([0,0.010,0.020])

plt.savefig(plotDir + 'SPE-Spectra.png',dpi=400)
plt.savefig(plotDir + 'SPE-Spectra.pdf',dpi=200)

plt.show()










