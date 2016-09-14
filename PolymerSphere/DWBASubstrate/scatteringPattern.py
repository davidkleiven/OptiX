import sys
sys.path.append("../../FresnelBEM")
sys.path.append("../../FresnelFDTD")
import matplotlib
import mplLaTeX as ml
matplotlib.rcParams.update( ml.params )
import numpy as np
from matplotlib import pyplot as plt
import json
import fresnelExact as fe

    
def dwba( waveVector, x, y, z, eps, mode="Born" ):
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
    firstDWBA = formFactorSphere(qR)*fe.rs(1.0,eps,1.0,1.0,angleWithSubstrateDeg, k)
    if ( mode == "First" ):
        return firstDWBA

    # Contrib from scattering after entering sphere
    qx = -kx_scat - waveVector[0]
    qR = np.sqrt( q_paralell**2 + qx**2 )
    angleWithSubstrateDeg = np.arccos(kx_scat/k)*180.0/np.pi
    secondDWBA = formFactorSphere(qR)*fe.rs(1.0,eps,1.0,1.0,angleWithSubstrateDeg,k)
    secondDWBA[-kx_scat>0.0] = 0.0 # Ray leaves sphere with kx = -kx_scat
    if ( mode == "Second" ):
        return secondDWBA

    # Contibution from scattering before enterin sphere and after
    qx = -kx_scat + waveVector[0]
    qR = np.sqrt( q_paralell**2 + qx**2 )
    angleWithSubstrateFirstDeg = np.arccos(np.abs(waveVector[0])/k)*180.0/np.pi
    angleWithSubstrateSecondDeg = np.arccos(kx_scat/k)*180.0/np.pi
    thirdDWBA = formFactorSphere(qR)*fe.rs(1.0,eps,1.0,1.0,angleWithSubstrateFirstDeg,k)*fe.rs(1.0,eps,1.0,1.0,angleWithSubstrateSecondDeg,k) 
    thirdDWBA[-kx_scat>0.0] = 0.0
    if ( mode == "Third" ):
        return thirdDWBA
    return formBorn + firstDWBA + secondDWBA + thirdDWBA
    
    
def formFactorSphere( qR ):
    form = ( np.sin(qR) - qR*np.cos(qR) )/qR**3
    return form

def grazingIncidence(grazingAngle, polarisation, epsSub):
    mu = 1.0
    epsInc = 1.0
    alpha_c = fe.criticalGrazingAngle( 1.0, np.sqrt(epsSub*mu) )
    print ("Critical angle: %.4f"%(alpha_c*180.0/np.pi))
    Rsphere = 1.0
    kR = 10.0
    k = np.zeros(3)
    k[0] = -kR*np.sin(grazingAngle*alpha_c)
    k[2] = kR*np.cos(grazingAngle*alpha_c)

    incAngle = 90.0 - grazingAngle*alpha_c*180.0/np.pi
    z = 100.0
    x = np.linspace(0.0,2.5*z,10001)
    born = dwba( k, x, 0.0, z, epsSub, mode="Born")
    f1 = dwba( k, x, 0.0, z, epsSub, mode="First")
    f2 = dwba( k, x, 0.0, z, epsSub, mode="Second")
    f3 = dwba( k, x, 0.0, z, epsSub, mode="Third")
    tot = dwba( k, x, 0.0, z, epsSub)
    alpha_f = np.arctan(x/z)/alpha_c

    colors = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99"]
    markers = ['-o', '-v', '-s', '-h', '-d']
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot( alpha_f, np.abs(born)**2, color=colors[0], label="A")
    ax.plot( alpha_f, np.abs(f1)**2, color=colors[1], label="B")
    ax.plot( alpha_f, np.abs(f2)**2, color=colors[2], label="C")
    ax.plot( alpha_f, np.abs(f3)**2, color=colors[3], label="D")
    ax.set_yscale('log')
    ax.set_xlabel("$\\alpha_f/\\alpha_c$")
    ax.set_ylabel("Intensity (a.u.)")
    if ( grazingAngle > 1.0 ):
        ax.text( 0.8*alpha_f[-1], 1E-10, "$\\alpha_i = %d\\alpha_c$"%(grazingAngle))
    elif( grazingAngle == 1 ):
        ax.text( 0.8*alpha_f[-1], 1E-10, "$\\alpha_i = \\alpha_c$") 
    else:
        ax.text( 0.8*alpha_f[-1], 1E-10, "$\\alpha_i = %.1f\\alpha_c$"%(grazingAngle))
    ax.legend( loc="lower right", frameon=False, ncol=4 )
    fname = "Figures/pattern1D.pdf"
    fig.savefig( fname, bbox_inches="tight" )
    print ("Figure written to %s"%(fname))

def main():
    grazingIncidence( 0.5, "TE", (1-1E-5+1j*1E-6)**2 )

if __name__ == "__main__":
    main() 
