# -*- coding: utf-8 -*-
import sys
import transmission as trans
import json
import h5py as h5
import numpy as np

def main( argv ):
    inpFile = ""
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            inpFile = arg.split("--file=")[1]
        elif ( arg.fin("--help") != -1 ):
            print ("Usage: python plotMultipleTransmission.py --file=<jsonfile>  [--help]")
            print ("file: Filename to the json file describing the datasets to plot")

    if ( inpFile == ""):
        print("No input file specified!")
        return

    try:
        infile = open(inpFile, 'r')
        entries = json.load(infile)
        infile.close()
    except Exception as exc:
        print ("Error when opening compare json file!")
        print str(exc)

    tPlot = trans.Transmission()
    tPlot.figname = entries["figname"]

    exitPlot = trans.ExitFields()
    exitPlot.figname = entries["exitfigname"]
    for entry in entries["entries"]:
        try:
            # Load json control file
            infile = open(entry["ctlfile"], 'r')
            stat = json.load(infile)
            infile.close()

            # Load data
            with h5.File(stat["Transmission"]["file"], 'r') as hf:
                data = np.array( hf.get( hf.keys()[0] ))
            z = np.linspace( stat["zDiscretization"]["min"], stat["zDiscretization"]["max"], len(data))

            if ( stat.get("waveguide").get("crd") ):
                crd = stat["waveguide"]["crd"]
            else:
                crd = "cartesian"

            if ( crd == "cylindrical" ):
                z *= stat["waveguide"]["RadiusOfCurvature"]
            tPlot.addData( z, data, entry["label"])

            # Load the far fields
            with h5.File(stat["farFieldFile"], 'r') as hf:
                dset = hf.get( "exitField" )
                xmin = dset.attrs.get("xmin")
                xmax = dset.attrs.get("xmax")
                data = np.array( dset )
            x = np.linspace(xmin, xmax, len(data))
            exitPlot.addData( x, data, entry["label"])
        except Exception as exc:
            print ("Error when reading data!")
            print str(exc)
    tPlot.plot()
    exitPlot.cutData()
    exitPlot.plot()

if __name__ == "__main__":
    main(sys.argv[1:])
