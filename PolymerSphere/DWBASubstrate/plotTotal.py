import sys
import grazing as graz
import numpy as np

MSG = "Usage: python plotTotal.py [--usefilm --help]+n"
MSG += "help - print this message\n"
MSG += "usefilm - if present the film reflection coefficients will be used\n"
MSG += "          otherwise the Fresnel coefficients will be used\n"
def main(argv):
    usefilm = False
    for arg in argv:
        if ( arg.find("--usefilm") != -1 ):
            usefilm = True
        elif ( arg.find("--help") != -1 ):
            print MSG
            return 0
    thickInUnitsOfR = 2.0
    angles = np.array([0.5, 1.0, 2.0,10.0])
    try:
        gz = graz.GrazingHandler(usefilm) 
        gz.setEpsilonSubst( (0.992+0.002j)**2 )
        gz.setFilmThickness(thickInUnitsOfR)
        gz.totalAngleSweep(angles)
    except Exception as e:
        print str(e)

if __name__ == "__main__":
    main(sys.argv[1:])
