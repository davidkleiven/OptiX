import sys
import grazingTransmission as gzt

MSG = "Usage: python plotTermsTransmission.py --d=<filmthick> --alpha=<incGrazAngle> [--help]\n"
MSG += "help: print this message\n"
MSG += "d: film thickness in units of the radius of the sphere\n"
MSG += "alpha: incident grazing angle in units of the critical angle\n"

def main(argv):
    dInUnitsOfR = 2.0
    alpha = 1.0
    n2 = 0.992+0.002j
    for arg in argv:
        if ( arg.find("--help") != -1 ):
            print MSG
            return 0
        elif ( arg.find("--d=") != -1 ):
            dInUnitsOfR = float(arg.split("--d=")[1])
        elif ( arg.find("--alpha=") != -1 ):
            alpha = float(arg.split("--alpha=")[1])
        else:
            print("Unknown argument %s"%(arg))
            return 0

    try:
        gh = gzt.GrazingTransmissionHandler()
        gh.setEpsilonSubst(n2**2)
        gh.setFilmThickness(dInUnitsOfR)
        gh.prepareDWBA(alpha)
        gh.plotTerms()
    except Exception as e:
        print str(e)
    return 0

if __name__ == "__main__":
    main(sys.argv[1:]) 
    
