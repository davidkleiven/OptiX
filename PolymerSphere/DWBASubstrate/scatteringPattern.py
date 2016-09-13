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

def formFactorSphere( incWaveVec, x, y, z, Rsphere ): 
    angle = np.arctan( np.sqrt(x**2 + y**2)/z )
    q = np.sqrt( np.sum(incWaveVec**2) )
    qR = q*Rsphere*np.sin(angle/2.0)
    qR[np.abs(qR) < 1E-5] = 1E-5
    form = ( np.sin(qR) - qR*np.cos(qR) )/qR**3
    return form/np.max(form)

def grazingIncidence(grazingAngle, polarisation, epsSub):
    mu = 1.0
    epsInc = 1.0
    alpha_c = fe.criticalGrazingAngle( 1.0, np.sqrt(epsSub*mu) )
    Rsphere = 1.0
    kR = 10.0
    k = np.zeros(3)
    k[1] = -kR*np.sin(grazingAngle*alpha_c)
    k[2] = kR*np.cos(grazingAngle*alpha_c)

    k_scattered = np.zeros(3)
    k_scattered[1] = kR*np.sin(grazingAngle*alpha_c)
    k_scattered[2] = kR*np.cos(grazingAngle*alpha_c)

    if ( polarisation == "TE" ):
        r = fe.rs( epsInc, epsSub, mu, mu, np.pi/2.0-grazingAngle*alpha_c, kR )
    elif ( polarisation == "TM" ):
        r = fe.rp( epsInc, epsSub, mu, mu, np.pi/2.0-grazingAngle*alpha_c, kR )
    else:
        print ("Unknown polarisation...")
        return
    
    z = 1000.0
    x = np.linspace( -z, z, 101 )
    y = np.linspace( -z, z, 101 )
    X, Y = np.meshgrid( x, y )
    f1 = formFactorSphere( k, X, Y, z, Rsphere )
    f2 = formFactorSphere( k_scattered, X, Y, z, Rsphere )
    f = np.abs(f1 + r*f2)**2
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.contourf(X,Y,f,200, map="gist_heat")
    
    fname = "Figures/pattern.png"
    fig.savefig( fname, bbox_inches="tight" )
    print ("Figure written to %s"%(fname))

    f1x = formFactorSphere( k, x, 0.0, z, Rsphere )
    f2x = formFactorSphere( k_scattered, x, 0.0, z, Rsphere )
    fx = np.abs(f1x + r*f2x)**2
    f1y = formFactorSphere( k, 0.0, y,z, Rsphere )
    f2y = formFactorSphere( k_scattered, 0.0, y, z, Rsphere )
    fy = np.abs(f1y + r*f2y)**2
     
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot( x, fx, 'k', label="Parallel" )
    ax.plot( y, fy, 'k--', label="Perpendicular")
    ax.set_xlabel("Distance from center of screen")
    ax.legend(frameon=False)
    fname = "Figures/pattern1D.pdf"
    fig.savefig( fname, bbox_inches="tight" )
    print ("Figure written to %s"%(fname))

def main():
    grazingIncidence( 0.0001, "TE", 1-1E-5+1j*1E-6 )

if __name__ == "__main__":
    main()
    
    
