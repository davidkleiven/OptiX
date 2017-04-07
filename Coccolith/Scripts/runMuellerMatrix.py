import sys
import subprocess
import json
import numpy as np
import h5py as h5
import multiprocessing as mp

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
    for 

def main():
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
