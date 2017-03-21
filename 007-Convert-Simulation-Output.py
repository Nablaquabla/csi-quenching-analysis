import h5py
import numpy as np

# ============================================================================
#  Convert LabView output into HDF5 format, split into signal and background
# ============================================================================
def convert_to_hdf5(mainDir,fileName):

    outFileName = '%s.h5'%(fileName.split('.')[0])
    # Create hdf5 file
    f = h5py.File(mainDir + outFileName,'w')

    # Add dataset / group to the hdf5 file
    keys = ['e_id','p_id','p_type','c_type','target','cell','dep_energy','time','x','y','z','weight','generation','n_scatter','inc_part','inc_energy']

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
    # Read csv file
    _tdata = np.loadtxt(mainDir + fileName).T
    cut = _tdata[0] <= 9e8
    dset = {}
    for i,k in enumerate(keys):
        dset[k] = f.create_dataset(k,data=_tdata[i][cut],dtype=datatypes[k],maxshape=(None,))
    del _tdata

    # Close the hdf5 file
    f.close()

# ============================================================================
#                                Run program
# ============================================================================
if __name__ == '__main__':
    mainDirectory = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Final-Data-Set/Simulations/Output/'
    angles = [45,39,33,27,24,21,18]
    for a in angles:
        convert_to_hdf5(mainDirectory,'%i.d'%a)



