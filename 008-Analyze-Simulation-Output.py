import numpy as np
import h5py
from scipy import interpolate
import matplotlib.pylab as plt

# Prepare quenching model functions for hydrogen and carbon used to convert recoils in EJ299
data_en_h, data_qf_h = np.loadtxt('./QF-H.dat',unpack=True,skiprows=1)
qf_h = interpolate.interp1d(data_en_h,data_qf_h/100.0,'linear',fill_value='extrapolate')

def qf_c(x):
    return 0.01

#xp = np.logspace(-2,4,1000)
#plt.plot(xp,qf_h(xp)*100.0,c='k')
#plt.show()


"""
      # What data type is necessary to save the information provided in the csv files. Mainly used to save some storage space
    datatypes = {'e_id': np.dtype(np.int32),
                 'p_id': np.dtype(np.int16),
                 'p_type': np.dtype(np.int16),
                 'c_type': np.dtype(np.int16),
                 'target': np.dtype(np.int32),
                 'cell': np.dtype(np.int16),
                 'dep_energy': np.dtype(np.float),
                 'time': np.dtype(np.float),
                 'x': np.dtype(np.float),
                 'y': np.dtype(np.float),
                 'z': np.dtype(np.float),
                 'weight': np.dtype(np.float),
                 'generation': np.dtype(np.int16),
                 'n_scatter': np.dtype(np.int16),
                 'inc_part': np.dtype(np.int16),
                 'inc_energy': np.dtype(np.float)}
"""

mainDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Final-Data-Set/Simulations/Output/'

def process(a):
    f = h5py.File(mainDir+'%i.h5'%a,'r')
    collType = f['c_type'][...]
    collCut = (collType == -99)

    eid = f['e_id'][...][collCut]
    target = f['target'][...][collCut]
    edep = f['dep_energy'][...][collCut]*1000.0
    cid = -1
    csiEnergies = []
    ejEnergy = 0
    ej_energies = np.array([])
    csi_energies = []
    for i in range(len(eid)):
        if eid[i] != cid:
            if i != 0:
                 ej_energies = np.concatenate((ej_energies,[ejEnergy]))
                 csi_energies.append(csiEnergies)
            csiEnergies = []
            ejEnergy = 0
            cid = eid[i]
        if target[i] == 1001:
            ejEnergy += qf_h(edep[i])*edep[i]
        elif target[i] == 6000:
            ejEnergy += qf_c(edep[i])*edep[i]
        else:
            csiEnergies.append(edep[i])
    f.close()

    noScatterings = np.zeros(20)
    for evt in csi_energies:
        if len(evt) < 20 and np.sum(evt) > 0:
            noScatterings[len(evt)] += 1

    maxLen = np.max(np.arange(20)[noScatterings > 0])
#    plt.step(np.arange(20),noScatterings,c='k')

    csi_edep = np.zeros((len(csi_energies),maxLen))
    for i in range(len(csi_energies)):
        for j in range(len(csi_energies[i])):
            csi_edep[i][j] = csi_energies[i][j]

    cut = np.sum(csi_edep,axis=1) > 0
    f = h5py.File(mainDir +'/%i-Sim-Results.h5'%a,'w')
    f.create_dataset('EJ-Edep', data=ej_energies[cut])
    f.create_dataset('CsI-Edep', data= csi_edep[cut,:])
    f.create_dataset('CsI-Scatters', data=noScatterings)
    f.close()

if __name__ == '__main__':
    for a in [45,39,33,27,24,21,18]:
        process(a)




















