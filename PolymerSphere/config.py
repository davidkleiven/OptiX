import sys

def main(argv):
    MSG = "Usage: python config.py --scuffh=<pathToScuffHeader> --buffh=<pathToBuffHeader> --libpath=<path to lib>\n"
    OUTFILE = "libHeaders.sh"

    outentries = {"scuffheader":"./", "buffheader":"./", "libpaths":"./"}
    for arg in argv:
        if ( arg.find("--help") != -1):
            print MSG
            return 0
        elif ( arg.find("--scuffh=") != -1 ):
            outentries["scuffheader"] = arg.split("--scuffh=")[1]
        elif ( arg.find("--buffh=") != -1 ):
            outentries["buffheader"] = arg.split("--buffh=")[1]
        elif ( arg.find("--libpath=") != -1 ):
            outentries["libpaths"] = arg.split("--libpath=")[1]
        else:
            print ("Unknown argument %s"%(arg))
            return 0

    ofile = open(OUTFILE, 'w')
    ofile.write("SCUFF_HEADER=%s\n"%(outentries["scuffheader"]))
    ofile.write("BUFF_HEADER=%s\n"%(outentries["buffheader"]))
    ofile.write("LIB_PATHS=%s\n"%(outentries["libpaths"]))
    ofile.close()
    return 0

if __name__ == "__main__":
    main(sys.argv[1:])
