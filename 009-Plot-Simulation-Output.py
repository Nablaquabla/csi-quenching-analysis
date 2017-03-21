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
#                          Custom plot functions
#=============================================================================#
def fill_between_step(x,y1,y2,color='r',alpha=0.25):
        if not len(x) == len(y1) == len(y2): return

        dx = np.diff(x)
        _tx = []
        _ty1 = []
        _ty2 = []

        for i, xx in enumerate(x):
            if i == 0:
                _tx.append(x[i])
                _ty1.append(y1[i])
                _ty2.append(y2[i])

                _tx.append(x[i]+dx[i]/2.0)
                _ty1.append(y1[i])
                _ty2.append(y2[i])

            elif (i == (len(x) - 1)):
                _tx.append(x[i]-dx[i-1]/2.0)
                _ty1.append(y1[i])
                _ty2.append(y2[i])

                _tx.append(x[i])
                _ty1.append(y1[i])
                _ty2.append(y2[i])

            else:
                _tx.append(x[i]-dx[i-1]/2.0)
                _ty1.append(y1[i])
                _ty2.append(y2[i])

                _tx.append(x[i+1]-dx[i]/2.0)
                _ty1.append(y1[i])
                _ty2.append(y2[i])

        plt.fill_between(_tx,_ty1,_ty2,color=color,alpha=alpha)
        return

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
mDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Final-Data-Set/Simulations/Output/'
plotDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Plots/'
saveFigs = True
for a in [18]:

    csiPlotMax = 30
    ej_en_th = 400

    f = h5py.File(mDir + '%i-Sim-Results.h5'%a,'r')

    ejedep = f['EJ-Edep'][...]/1000.0
    csiedep = f['CsI-Edep'][...]
    csiscatters = f['CsI-Scatters'][...]

    plt.figure(figsize=(6.5,4.0),edgecolor='k',facecolor='w')
    gs = gridspec.GridSpec(2,2,height_ratios=[1,1.5],width_ratios=[3,1])
    gs.update(left=0.10, right=0.9, top=0.9, bottom=0.13, wspace=0.02, hspace=0.075)

    ax = []

    ax.append(plt.subplot(gs[0,0]))
    ax.append(plt.subplot(gs[1,0]))
    ax.append(plt.subplot(gs[0,1]))
    ax.append(plt.subplot(gs[1,1]))

    # Upper left
    plt.sca(ax[0])
    n,b = np.histogram(ejedep,201,[-0.005,2.005])
    x = b[:-1] + 0.5*np.diff(b)
    plt.step(x,n,where='mid',c=colors.black)
#    plt.xlabel(r'EJ299 Energy [keV$_\mathrm{ee}$]')
    plt.ylabel(r'Counts / 10 keV$_\mathrm{ee}$',labelpad=-3)
    plt.yscale('log',nonposy='clip')
    plt.tick_params(axis='x',labelbottom='off')
    plt.xlim(0,1.750)
    plt.xticks([0,0.500,1.000,1.500])
    plt.ylim(5,5e3)

    # Lower left
    plt.sca(ax[1])
    x = ejedep
    y = np.sum(csiedep,axis=1)
    plt.scatter(x,y,marker='+',alpha=0.2,color=colors.black)
    plt.axvline(ej_en_th/1000.0,color=colors.red,linestyle='dashed')
    plt.xlim(0,1.750)
    plt.xticks([0,0.500,1.000,1.500])
    plt.ylim(0,csiPlotMax)
    plt.xlabel(r'EJ299 Energy [MeV$_\mathrm{ee}$]')
    plt.ylabel(r'E$_\mathrm{CsI}$ [keV$_\mathrm{nr}$]',labelpad=6)

    # Upper right
    plt.sca(ax[2])
    b = (np.arange(20) > -1) * (np.arange(20) < 6)
    plt.errorbar(np.arange(20)[b],csiscatters[b],yerr=np.sqrt(csiscatters[b]),color=colors.black,linestyle='None',marker='.')
    plt.step(np.arange(20)[b],csiscatters[b],color=colors.black,where='mid')
    print f['CsI-Scatters'][1]
    print np.sum(f['CsI-Scatters'][2:5])
    print f['CsI-Scatters'][1]/np.sum(f['CsI-Scatters'][...])
    plt.xlim(0.5,4.5)
    plt.ylim(1,1e6)
    plt.yscale('log',nonposy='clip')
    plt.ylabel('# of Histories',rotation=270,labelpad=16)
    plt.xlabel('Scatters in CsI[Na]',labelpad=13)
    plt.yticks([10,1e3,1e5])
    plt.xticks([1,2,3,4])
    plt.tick_params(axis='y',labelleft='off',labelright='on')
    plt.tick_params(axis='x', which='major', pad=-20,labelbottom='off',labeltop='on')
    ax[2].yaxis.set_label_position("right")
    ax[2].xaxis.set_label_position("top")

    # Lower right
    plt.sca(ax[3])
    eje = f['EJ-Edep'][...]
    y = np.sum(f['CsI-Edep'],axis=1)

    c = eje > ej_en_th
    n,b = np.histogram(y[c],121,[-0.125,30.125])
    x = b[:-1] + 0.5*np.diff(b)
    plt.step(n,x,where='mid',c=colors.red,label='Cut')
    x2,pars,xfit,yfit = ef.fit('gauss',x,n,'None',[np.max(n),x[np.argmax(n)],0.25,0])
    plt.plot(yfit,xfit,c=colors.blue,ls='dashed')
    print pars[0]
    n,b = np.histogram(y,121,[-0.125,30.125])
    x = b[:-1] + 0.5*np.diff(b)
    plt.step(n,x,where='mid',c=colors.black,label='Full')
    plt.legend(loc='best',fontsize=10,framealpha=0)
    plt.xticks([0,400,800,1200,1600])
    plt.ylabel(r'E$_\mathrm{CsI}$ [keV$_\mathrm{nr}$]',rotation=270,labelpad=24)
    plt.xlabel('Cts / 0.25 keV$_\mathrm{nr}$')
    plt.tick_params(axis='y',labelleft='off',labelright='on')
    ax[3].yaxis.set_label_position("right")
    plt.ylim(0,csiPlotMax)
    if saveFigs:
        plt.savefig(plotDir + '%iDeg-Sim-Out.png'%a,dpi=400)
        plt.savefig(plotDir + '%iDeg-Sim-Out.pdf'%a,dpi=200)
    f.close()
    plt.show()


















