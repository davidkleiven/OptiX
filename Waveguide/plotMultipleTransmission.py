import sys
import transmission as trans
import json
import h5py as h5

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
        print str(exc)

    tPlot = trans.Transmission()
    for entry in entries["entries"]:
        try:
            # Load json control file
            infile = open(entries["ctlfile"], 'r')
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
            tPlot.addSet( z, data, entry["label"])
        except Exception as exc:
            print str(exc)
    tPlot.plot()

if __name__ == "__main__":
    main(sys.argv[1:])
