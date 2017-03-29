import sys
import numpy as np
import multiprocessing as mp
import subprocess
import json
import os
import time
import h5py as h5

TEMP_DIR = "data/temporaryBandStructure"

class BandDiagramComputer:
    def __init__( self ):
        self.jsontemplate = None
        self.kx = 0.0
        self.ky = 0.0
        self.uid = 0
        self.blochPathIndx = 0

def computeBandDiagram( bandDiagComp ):
    # Create input file for the simulation
    bandDiagComp.jsontemplate["bloch"] = [bandDiagComp.kx,bandDiagComp.ky]
    bandDiagComp.jsontemplate["prefix"] = TEMP_DIR + "/bandStrucPythonParallel_bpath%d"%(bandDiagComp.blochPathIndx)
    bandDiagComp.jsontemplate["uid"] = bandDiagComp.uid
    outname = TEMP_DIR+"/inputfile%d.json"%(bandDiagComp.uid)
    ofile = open(outname,'w')
    json.dump( bandDiagComp.jsontemplate, ofile )
    ofile.close()

    logfile = open(TEMP_DIR+"/outfile%d.log"%(bandDiagComp.uid), 'w')
    # Launch the C++ executable with the newly created input file
    subprocess.call( ["./bandStructure.out", outname], stdout=logfile )
    logfile.close()
    print ("Job %d finished"%(bandDiagComp.uid))

def collectHDF5():
    isFirst = True
    tstamp = time.strftime("%Y%m%d_%H%M")
    hfile = h5.File( "data/bandStructure_%s.h5"%(tstamp) )
    counter = 0
    for fn in os.listdir(TEMP_DIR):
        if ( fn.find(".h5") == -1 ):
            continue

        with h5.File(TEMP_DIR+"/"+fn, 'r') as hf:
            if ( isFirst ):
                # Store some extra datasets here
                eps = np.array( hf.get("eps") )
                hfile.create_dataset("eps", data=eps)
                specFreqs = np.array( hf.get("freq") )
                hfile.create_dataset("freq", data=specFreqs )
                domain = np.array( hf.get("domain") )
                hfile.create_dataset("domain", data=domain)
                srcPos = np.array( hf.get("sourcePos") )
                hfile.create_dataset("sourcePos", data=srcPos)
                isFirst = False
            group = hfile.create_group("Run%d"%(counter))
            counter += 1
            bloch = np.array( hf.get("bloch") )
            group.attrs["kx"] = bloch[0]
            group.attrs["ky"] = bloch[1]
            ldos = np.array( hf.get("ldos") )
            group.create_dataset( "ldos", data=ldos)
            ez = np.array( hf.get("Ez") )
            group.create_dataset("Ez", data=ez)
            bpathPos = fn.find("bpath")
            pathNum = int(fn[bpathPos+5])
            group.attrs["blochPath"] = pathNum
    hfile.close()


def main( argv ):
    if ( len(argv) != 2 ):
        print ("Usage: python launchBandStructureParalle.py inpufile.json nprocessors")
        return
    subprocess.call( ["mkdir", "-p", TEMP_DIR] )
    nproc = int(argv[1])
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
            newRun.blochPathIndx = i
            uid += 1
            bloch.append(newRun)

    print ("Total number of runs %d"%(len(bloch)))
    workers = mp.Pool(nproc)
    workers.map( computeBandDiagram, bloch )
    collectHDF5()

    # Clean up
    subprocess.call(["rm", "-r", TEMP_DIR])

if __name__ == "__main__":
    main( sys.argv[1:] )
