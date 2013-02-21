from synapse_simulation import synapse_simulation
import numpy as np
import matplotlib.pyplot as plot
from mol_type import mol_type
from simulation_values import simulation_values
from scipy.optimize.optimize import brute
from scipy.optimize import fmin
from sys import stdout
from datetime import datetime
from norm import norm
import os
import time
import shutil
import glob

def run_simulation(speeds):
    a = synapse_simulation(speeds)
    myosins = np.zeros(simulation_values.num_steps)
    actins = np.zeros(simulation_values.num_steps)
    a.initialize()

    # Declares variable that tracks the clocktime between steps
    elapsed = float(0)
    # Helps make sure pictures are named with proper zeros
    digits = len(str(simulation_values.num_steps-1))
    namedir = make_directory()
    
    for i in range(simulation_values.num_steps):

        if (simulation_values.ok_count > 1):
            stdout.write("\r%d%% complete\t%s molecules\t%f sec\t%f band" % \
                         ((int(float(i)/float(simulation_values.num_steps)*100)),\
                          a.how_many,elapsed,a.band))
        else:
            stdout.write("\r%d%% complete\t%s molecules\t%f sec\t%s band" % \
                         ((int(float(i)/float(simulation_values.num_steps)*100)),\
                          a.how_many,elapsed,a.band))
        stdout.flush()

        start = time.time()

        a.step()
        #a.dist_to_band()            
        myosins[i] = a.myo_count
        actins[i] = a.act_count

        # Generates pictures
        mmyo,mact,mok = [0,0,0]
        plot.clf()
        for molecule in a.mol_all:
            c = "g"
            m = "o"
            size=1000./float(simulation_values.carrying_capacity)
            if (molecule.t == mol_type.ACTIN):
                if mact == 0:
                    lab = "Actin"
                else:
                    lab = "_nolegend_"
                mact += 1
            if (molecule.t == mol_type.OKT3):
                c = "k"
                size = 200.
                if mok == 0:    
                    lab = "OKT3"
                else:
                    lab = "_nolegend_"
                mok += 1
            if (molecule.t == mol_type.MYOSIN):
                c = "r"
                if mmyo == 0:
                    lab = "Myosin"
                else:
                    lab = "_nolegend_"
                mmyo += 1
            plot.scatter([molecule.pos[0]],[molecule.pos[1]],s=size,color=c,marker=m,label=lab)
        plot.axis([0,100,0,100])
        plot.legend(loc="upper left")
        if (len(str(i)) < digits):
            picname = "0%s" % i
        else:
            picname = "%s" % i
        plot.savefig(str(namedir+picname))
        
        # Plots and saves population curves    
        if (i == simulation_values.num_steps-1):
            plot.figure()
            plot.axis([0, (simulation_values.num_steps*simulation_values.time_step), 0, 2./5.*simulation_values.carrying_capacity])
            plot.plot(np.arange(simulation_values.num_steps)*simulation_values.time_step,myosins,color='r',label="Myosin")
            plot.plot(np.arange(simulation_values.num_steps)*simulation_values.time_step,actins,color='g',label="Actin")
            plot.legend(loc="upper left")
            plot.grid(True)
            plot.savefig(namedir+"populations")

        elapsed = time.time() - start

    print
    '''
    hist_mol = np.array([])
    hist_mol = a.mol_all[4].hist.transpose()
    dist_band = np.abs(50 - hist_mol[1])
    veloc = np.diff(hist_mol)
    veloc = veloc.transpose()
    speeds = np.zeros(veloc.shape[0])
    for i in range(veloc.shape[0]):
        speeds[i] = norm(veloc[i])
    speeddist = np.array([dist_band[0:(dist_band.size-1)],speeds])
    speeddist =  speeddist.transpose()
    #speeddist = np.sort(speeddist.view('i8,i8'), order=['f0'], axis=0).view(int)
    plot.figure()
    plot.scatter(dist_band[0:(dist_band.size-1)],speeds)
    #plot.figure()
    #plot.scatter(speeddist[0],speeddist[1])
    plot.savefig(namedir+"speed")
    '''

def make_directory():
    """Makes directory where pictures are saved and puts all the .py
    files in that directory
    """

    def make_dirname():
        now = datetime.today()    
        return "pictures/%sm%sd%sy%s%s/" % \
                  (str(now.month),str(now.day),str(now.year),str(now.hour),str(now.minute))

    dirname = make_dirname()

    while os.path.exists(dirname):
        print "Directory exists. Waiting 5 seconds to retry..."
        time.sleep(5)
        dirname = make_dirname()

    os.makedirs(dirname)

    # Copies in the .py files for reference
    py_files = glob.glob("./*.py")
    for i in py_files:
        shutil.copy(i,dirname)
    
    return dirname


run_simulation([])
