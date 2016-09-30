import sys
import grazing as graz
import grazingTransmission as grazT
import numpy as np
import transformDetector as td

MSG = "Usage: python plotTotal.py [--usefilm --help --transmission]\n"
MSG += "help - print this message\n"
MSG += "usefilm - if present the film reflection coefficients will be used\n"
MSG += "          otherwise the Fresnel coefficients will be used\n"
MSG += "transmission - if present the field transmitted through the film will be used\n"
def main(argv):
    usefilm = False
    useTransmission = False
    for arg in argv:
        if ( arg.find("--usefilm") != -1 ):
            usefilm = True
        elif ( arg.find("--help") != -1 ):
            print MSG
            return 0
        elif ( arg.find("--transmission") != -1 ):
            useTransmission = True
    thickInUnitsOfR = 2.0
    angles = np.array([0.5, 1.0, 2.0,5.0])
    try:
        if ( useTransmission ):
            gz = grazT.GrazingTransmissionHandler()
        else:
            gz = graz.GrazingHandler(usefilm)
        gz.setEpsilonSubst( (0.992+0.002j)**2 )
        gz.setFilmThickness(thickInUnitsOfR)
        gz.detectorTransform = td.DetectorCenterBeam() # Plot in terms of deflection angle
        if ( useTransmission ):
            gz.plotPattern(angles)
        else:
            gz.totalAngleSweep(angles)
    except Exception as e:
        print str(e)

if __name__ == "__main__":
    main(sys.argv[1:])
