import os
import multiprocessing as mp
import time

#outputDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Final-Data-Set/Simulations/Output/'
#inputDir = 'F:/Work-Data-Storage/CsI/CsI-Quenching/Final-Data-Set/Simulations/Input/'
outputDir = './Output/'
inputDir = './Input/'

def run(a):

    # Delete old output files
    if os.path.exists('%s%i.o'%(outputDir,a)):
        del_o = 'rm %s%i.o'%(outputDir,a)
        os.system(del_o)
    if os.path.exists('%s%i.d'%(outputDir,a)):
        del_d = 'rm %s%i.d'%(outputDir,a)
        os.system(del_d)

    # Prepare command
    cmd = 'MCNPX-Polimi-2.exe i=%s%i.i o=%s%i.o d=%s%i.d'%(inputDir,a,outputDir,a,outputDir,a)
#    print cmd
    os.system(cmd)

def update_input(a):
    # Define rotation cards for all angles
    tr_cards = {18: 'tr1 104.509 34.3204 0 0.950081 0.312003 0 -0.312003 0.950081 0 0 0 1',
                21: 'tr1 92.9003 37.0071 0 0.929003 0.370071 0 -0.370071 0.929003 0 0 0 1',
                24: 'tr1 91.5663 40.1948 0 0.915663 0.401948 0 -0.401948 0.915663 0 0 0 1',
                27: 'tr1 89.1323 45.3368 0 0.891323 0.453368 0 -0.453368 0.891323 0 0 0 1',
                33: 'tr1 75.5231 48.9516 0 0.839146 0.543907 0 -0.543907 0.839146 0 0 0 1',
                39: 'tr1 69.9431 56.6388 0 0.777146 0.62932 0 -0.62932 0.777146 0 0 0 1',
                45: 'tr1 56.8148 56.3212 0 0.710185 0.704015 0 -0.704015 0.710185 0 0 0 1'}

    # Read source definition from template file
    sdef = []
    with open('sim-sdef.txt','r') as fSource:
        for line in fSource:
            sdef.append(line)

    with open('%s%i.i'%(inputDir,a),'w') as fOut:
        with open('./sim-duke-template.i','r') as fIn:
            for line in fIn:
                if 'PYTHON-REPLACE-TRCARD' in line:
                    fOut.write(tr_cards[a])
                    fOut.write('\n')
                elif 'PYTHON-REPLACE-SDEF' in line:
                    for SL in sdef:
                        fOut.write(SL)
                    fOut.write('\n')
                else:
                    fOut.write(line)


if __name__ == '__main__':

#    angles = [45,39]
#    angles = [33,27,24]
#    angles = [21,18]
    angles = [45,39,33,18,21,24,27]
#    angles = [45,39,33]

    procs = []
    for a in angles:
        update_input(a)
        p = mp.Process(target=run, args=(a,))
        p.start()
        procs.append(p)
        time.sleep(2)