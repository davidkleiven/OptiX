import sys
sys.path.append("../../FresnelBEM")
sys.path.append("../../FresnelFDTD")
sys.path.append("../../")
import matplotlib
import mplLaTeX as ml
matplotlib.rcParams.update( ml.params )
import numpy as np
from matplotlib import pyplot as plt
import json
import fresnelExact as fe
import colorScheme as cs
import scatteringStructures as scat
import grazing as graz

    
HELP_MSG = "Usage: python scatteringPattern.py [--usefilm --d=filmthickness --alpha=0.5]\n"
HELP_MSG += "--help - print this message\n"
HELP_MSG += "--usefilm - use the thinfilm reflection coefficients in stead\n"
HELP_MSG += "--d - film thickness in units of the radius of the sphere\n"
HELP_MSG += "--alpha - incident angle in units of the critical angle, 1 is the default\n"

def dwba( waveVector, x, y, z, scatObj, mode="Born" ):
    k = np.sqrt( np.sum(waveVector**2) )
    q_paralell = np.sqrt( waveVector[1]**2 + waveVector[2]**2 ) 
    rHat_x = x/np.sqrt(x**2+y**2+z**2)
    kx_scat = k*rHat_x
    qx = kx_scat - waveVector[0]
    qR = np.sqrt( q_paralell**2 + qx**2 )
    qR[qx<0.0] = -qR[qx<0.0]
    formBorn = formFactorSphere(qR)
    if ( mode == "Born" ):
        return formBorn
    
    # Contrib from scattered before entering sphere
    qx = kx_scat + waveVector[0]
    qR = np.sqrt( q_paralell**2 + qx**2 )
    angleWithSubstrateDeg = np.arccos(-waveVector[0]/k)*180.0/np.pi
    firstDWBA = formFactorSphere(qR)*scatObj.reflection(angleWithSubstrateDeg)#fe.rs(1.0,eps,1.0,1.0,angleWithSubstrateDeg, k)
    if ( mode == "First" ):
        return firstDWBA

    # Contrib from scattering after entering sphere
    qx = -kx_scat - waveVector[0]
    qR = np.sqrt( q_paralell**2 + qx**2 )
    angleWithSubstrateDeg = np.arccos(kx_scat/k)*180.0/np.pi
    secondDWBA = formFactorSphere(qR)*scatObj.reflection(angleWithSubstrateDeg)#fe.rs(1.0,eps,1.0,1.0,angleWithSubstrateDeg,k)
    secondDWBA[-kx_scat>0.0] = 0.0 # Ray leaves sphere with kx = -kx_scat
    if ( mode == "Second" ):
        return secondDWBA

    # Contibution from scattering before enterin sphere and after
    qx = -kx_scat + waveVector[0]
    qR = np.sqrt( q_paralell**2 + qx**2 )
    angleWithSubstrateFirstDeg = np.arccos(np.abs(waveVector[0])/k)*180.0/np.pi
    angleWithSubstrateSecondDeg = np.arccos(kx_scat/k)*180.0/np.pi
    thirdDWBA = formFactorSphere(qR)*scatObj.reflection(angleWithSubstrateFirstDeg)*scatObj.reflection(angleWithSubstrateSecondDeg)#fe.rs(1.0,eps,1.0,1.0,angleWithSubstrateFirstDeg,k)*fe.rs(1.0,eps,1.0,1.0,angleWithSubstrateSecondDeg,k) 
    thirdDWBA[-kx_scat>0.0] = 0.0
    if ( mode == "Third" ):
        return thirdDWBA
    return formBorn + firstDWBA + secondDWBA + thirdDWBA
    
    
def formFactorSphere( qR ):
    form = ( np.sin(qR) - qR*np.cos(qR) )/qR**3
    return form

def grazingIncidence(grazingAngle, polarisation, epsSub, useFilm=False, dInUnitsOfR=None):
    mu = 1.0
    epsInc = 1.0
    if ( useFilm ):
        scatObj = scat.FilmScatterer()
        if ( dInUnitsOfR is None ):
            print ("Cannot use film coefficients without a thickness")
            return
        scatObj.thickness = dInUnitsOfR 
    else:
        scatObj = scat.PlaneScatterer()

    scatObj.eps2 = epsSub
    alpha_c = fe.criticalGrazingAngle( 1.0, np.sqrt(epsSub*mu) )
    print ("Critical angle: %.4f"%(alpha_c*180.0/np.pi))
    Rsphere = 1.0
    kR = 5.0
    scatObj.k = kR
    k = np.zeros(3)
    k[0] = -kR*np.sin(grazingAngle*alpha_c)
    k[2] = kR*np.cos(grazingAngle*alpha_c)

    incAngle = 90.0 - grazingAngle*alpha_c*180.0/np.pi
    z = 100.0
    x = np.linspace(0.0,2.5*z,10001)
    born = dwba( k, x, 0.0, z, scatObj, mode="Born")
    f1 = dwba( k, x, 0.0, z, scatObj, mode="First")
    f2 = dwba( k, x, 0.0, z, scatObj, mode="Second")
    f3 = dwba( k, x, 0.0, z, scatObj, mode="Third")
    tot = dwba( k, x, 0.0, z, scatObj)
    alpha_f = np.arctan(x/z)/alpha_c

    markers = ['-o', '-v', '-s', '-h', '-d']
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot( alpha_f, np.abs(born)**2, color=cs.COLORS[0], label="A")
    ax.plot( alpha_f, np.abs(f1)**2, color=cs.COLORS[1], label="B")
    ax.plot( alpha_f, np.abs(f2)**2, color=cs.COLORS[2], label="C")
    ax.plot( alpha_f, np.abs(f3)**2, color=cs.COLORS[3], label="D")
    ax.set_yscale('log')
    ax.set_xlabel("$\\alpha_f/\\alpha_c$")
    ax.set_ylabel("Intensity (a.u.)")
    if ( grazingAngle > 1.0 ):
        ax.text( 0.8*alpha_f[-1], 1E-5, "$\\alpha_i = %d\\alpha_c$"%(grazingAngle))
    elif( grazingAngle == 1 ):
        ax.text( 0.8*alpha_f[-1], 1E-5, "$\\alpha_i = \\alpha_c$") 
    else:
        ax.text( 0.8*alpha_f[-1], 1E-5, "$\\alpha_i = %.1f\\alpha_c$"%(grazingAngle))
    ax.legend( loc="lower right", frameon=False, ncol=4 )
    
    if ( useFilm ):
        fname = "Figures/patternFilm1D.pdf"
    else:
        fname = "Figures/pattern1D.pdf"

    fig.savefig( fname, bbox_inches="tight" )
    print ("Figure written to %s"%(fname))

def main(argv):
    useFilm = False
    dInUnitsOfR = None
    alpha = 1.0
    for arg in argv:
        if ( arg.find("--usefilm") != -1 ):
            useFilm = True
        elif ( arg.find("--d=") != -1 ):
            dInUnitsOfR = float(arg.split("--d=")[1])
        elif ( arg.find("--help") != -1 ):
            print HELP_MSG
            return 0
        elif ( arg.find("--alpha=") != -1 ):
            alpha = float( arg.split("--alpha=")[1])
    #grazingIncidence( alpha, "TE", (0.992+1j*0.002)**2, useFilm=useFilm, dInUnitsOfR=dInUnitsOfR )

    # Test the new GrazingHandler class
    gz = graz.GrazingHandler(useFilm)
    gz.setEpsilonSubst( (0.992+1j*0.002)**2 )
    gz.setFilmThickness(dInUnitsOfR)
    gz.prepareDWBA(alpha)
    gz.plotTerms()

if __name__ == "__main__":
    main(sys.argv[1:]) 
