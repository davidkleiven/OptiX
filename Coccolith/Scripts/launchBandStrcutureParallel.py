import sys
import numpy as np
import multiprocessing as mp
import subprocess
import json

TEMP_DIR = "data/temporaryBandStructure"

class BandDiagramComputer:
    def __init__( self ):
        self.jsontemplate = None
        self.kx = 0.0
        self.ky = 0.0
        self.uid = 0

def computeBandDiagram( bandDiagComp ):
    # Create input file for the simulation
    bandDiagComp.jsontemplate["bloch"] = [bandDiagComp.kx,bandDiagComp.ky]
    bandDiagComp.jsontemplate["prefix"] = TEMP_DIR + "/bandStrucPythonParallel"
    bandDiagComp.jsontemplate["uid"] = bandDiagComp.uid
    outname = TEMP_DIR+"/inputfile%d.json"%(bandDiagComp.uid)
    ofile = open(outname,'w')
    json.dump( bandDiagComp.jsontemplate, ofile )
    ofile.close()

    # Launch the C++ executable with the newly created input file
    subprocess.call( ["./bandStructure.out", outname] )

def main( argv ):
    subprocess.call( ["mkdir", "-p", TEMP_DIR] )

    # Parse the inputfile
    infile = open( argv[0], 'r' )
    params = json.load( infile )
    infile.close()

    infile = open( params["inputfileTemplate"], 'r' )
    templateInput = json.load(infile)
    infile.close()

    # Create list of pairs of Bloch-vectors
    bloch = []
    uid = 0
    for i in range(0, len(params["blochvectorPath"])-1):
        kx = params["blochvectorPath"][i]["kx"]
        ky = params["blochvectorPath"][i]["ky"]
        kxEnd = params["blochvectorPath"][i+1]["kx"]
        kyEnd = params["blochvectorPath"][i+1]["ky"]
        dkx = (kxEnd-kx)/params["numberOfKVecInEachInterval"]
        dky = (kyEnd-ky)/params["numberOfKVecInEachInterval"]
        for j in range(0, params["numberOfKVecInEachInterval"]):
            newRun = BandDiagramComputer()
            newRun.jsontemplate = templateInput
            newRun.kx = kx+j*dkx
            newRun.ky = ky+j*dky
            newRun.uid = uid
            uid += 1
            bloch.append(newRun)

    workers = mp.Pool(4)
    workers.map( computeBandDiagram, bloch )

    # Clean up
    subprocess.call(["rm", "-r", TEMP_DIR])

if __name__ == "__main__":
    main( sys.argv[1:] )
