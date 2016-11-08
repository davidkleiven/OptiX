import numpy as np
import sys

from scipy import stats
from scipy import signal as sig
import json
import angleComputer as ac
import wpdExtract as wpd

def main( argv ):
    fname = ""
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        elif ( arg.find("--help") != -1 ):
            print ("Usage: python computeAngles.py --file=<json-file>")
            return 1

    if ( fname == "" ):
        print ("No json file specified!")
        return 1

    try:
        infile = open(fname, 'r')
        data = json.load(infile)
        infile.close()
    except Exception as exc:
        print ("Error when loading json file!")
        print (str(exc))
        return 1

    try:
        wpdEx = wpd.WebPlotDigitizerDsetExtractor()
        wpdEx.parse( data )
    except Exception as exc:
        print ("Problem during wpd parsing")
        print str(exc)
        return 1

    angC = ac.AngleComputer()
    angC.wpd = wpdEx
    angC.plotView()

if __name__ == "__main__":
    main(sys.argv[1:])
