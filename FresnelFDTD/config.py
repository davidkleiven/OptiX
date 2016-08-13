# This is a configuration file that sets the MEEP path
import sys

OUTFILE = "paths.sh"
def main(argv):
    if ( len(argv) != 2 ):
        print ("Usage: python config.py --meepdir=/path/to/meep/src --meeplibdir=/path/to/meep/lib")
        return 1

    meepdir = ""
    meeplibdir = ""
    for arg in argv:
        if ( arg.find("--meepdir=") != -1 ):
            meepdir = arg.split("--meepdir=")[1]
        elif ( arg.find("--meeplibdir=") != -1 ):
            meeplibdir = arg.split("--meeplibdir=")[1]

    if ( meepdir == "" ):
        print ("Did not find meepdir in the arguments...")
        return 1

    if ( meeplibdir == "" ):
        print ("Did not find meeplibdir in the arguments...")
        return 1

    out = open( OUTFILE, 'w' )
    out.write("MEEP_IDIR=\"%s\"\n"%(meepdir)) 
    out.write("MEEP_LDIR=\"%s\"\n"%(meeplibdir))
    out.close()
    return 0

if __name__ == "__main__":
    main(sys.argv[1:])
