import sys
import subprocess
import json
import numpy as np
import h5py as h5
import multiprocessing as mp
import time

tempDir = "data/muellerCalculation"
templateInputFile = "InputFiles/visibleLightCaCO3ne.json"

class MuellerMatrixSimulation:
    def __init__( self ):
        self.uid = 0
        self.stokesVector = []

def launch( muellerSim ):
    # Load input file
    inf = open( templateInputFile, 'r' )
    tmpInput = json.load(inf)
    inf.close()

    tmpInput["incStokesVector"] = muellserSim.stokesVector
    tmpInput["prefix"] = "muellerCalculation/CaCO3_90deg"
    ofname = tempDir+"/input%d.json"%(muellerSim.uid)
    out = open(ofname, 'w')
    json.dump(tmpInput,out)
    out.close()

    subprocess.call(["mpiexec", "-n", "12", "./coccolith.out", ofname])

def mergeHDF5s():
    tstamp = time.strftime("%Y%m%d_%H%M")
    hfile = h5.File( "data/muellerMatrix_%s.h5"%(tstamp) )
    isFirst = True
    counter = 0
    for fname in os.listdir(tempDir):
        if (fname.find(".h5") == -1 ):
            continue

        with h5.File(tempDir+"/"+fname) as hf:
            if ( isFirst ):
                hfile.create_dataset("eps", data=np.array(hf.get("eps")))
                isFirst = False

            group = hfile.create_group("Run%d"%(counter))
            counter += 1
            group.create_dataset("freq", data=np.array( hf.get("freq") ) )

            # Collect everything that contains stokes or Stokes
            for key in hf.keys():
                if (( "stokes" in key ) or ("Stokes" in key) or ("stoke" in key) or ("Stoke" in key)):
                    group.create_dataset(key,data=np.array(hf.get(key)))
        # Clean up
        subprocess.call(["rm", tempDir+"/"+fname])
    hfile.close()

def main( argv ):
    global templateInputFile
    if  ( len(argv) != 1 ):
        print ("Usage: python InputFileTemplate.json")
        return
    templateInputFile = argv[0]
    
    # Create a temporary directory
    subprocess.call(["mkdir", "-p", tempDir])

    stokesVectors=[[1,1,0,0],[1,-1,0,0],[1,0,0,1],[1,0,1,0]]

    stokesSimulations = []
    for i in range(0,len(stokesVectors)):
        newstoke = MuellerMatrixSimulation()
        newstoke.uid = i
        newstoke.stokesVector = stokesVectors[i]

    workers = mp.Pool(2)
    workers.map( launch, stokesSimulations )
    mergeHDF5s()
