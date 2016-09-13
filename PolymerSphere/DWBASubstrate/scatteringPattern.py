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

def Born1( incWaveVec, x, y, z, Rsphere ):
    q_parallel = qPar( incWaveVec )
    k_scattered_z = kScatteredPerp( incWaveVec, x, y, z )
    q_norm = k_scattered_z - incWaveVec[2]
    return formFactorSphere( q_parallel, q_norm, Rsphere )
 
def reflectedFromSubsrate( incWaveVec, x, y, z, Rsphere, epsSubst, pol="TE" ):
    q_parallel = qPar( incWaveVec )
    k_scattered_z = kScatteredPerp( incWaveVec, x, y, z )
    q_norm = k_scattered_z + incWaveVec[2]
    form = formFactorSphere( q_parallel, q_norm, Rsphere )
    k = np.sqrt( np.sum(incWaveVec**2) )
    if ( pol == "TE" ):
        return form*fe.rs( 1.0, epsSubst, 1.0, 1.0, incAngle(incWaveVec), k )
    return form*fe.rp( 1.0, epsSubst, 1.0, 1.0, incAngle(incWaveVec), k )

def reflectedFromSubstrateAfter( incWaveVec, x, y, z, Rsphere, epsSubst, pol="TE" ):
    q_parallel = qPar( incWaveVec )
    k_scattered_z = kScatteredPerp( incWaveVec, x, y, z )
    q_norm = -k_scattered_z - incWaveVec[2]
    form = formFactorSphere( q_parallel, q_norm, Rsphere )
    k = np.sqrt( np.sum(incWaveVec**2) )
    angle = np.arccos(k_scattered_z/k)*180.0/np.pi
    if ( pol == "TE" ):
        return form*fe.rs( 1.0, epsSubst, 1.0, 1.0, angle, k )
    return form*fe.rp( 1.0, epsSubst, 1.0, 1.0, angle, k )
    
def reflectedFromSubsrateBeforeAndAfter( incWaveVec, x, y, z, Rsphere, epsSubst, pol="TE" ):
    q_parallel = qPar( incWaveVec )
    k_scattered_z = kScatteredPerp( incWaveVec, x, y, z )
    q_norm = -k_scattered_z + incWaveVec[2]
    form = formFactorSphere( q_parallel, q_norm, Rsphere )
    k = np.sqrt( np.sum(incWaveVec**2) )
    angle = np.arccos(k_scattered_z/k)*180.0/np.pi
    if ( pol == "TE" ):
        return form*fe.rs( 1.0, epsSubst, 1.0, 1.0, incAngle(incWaveVec), k )*fe.rs(1.0, epsSubst, 1.0,1.0, angle, k)
    return form*fe.rp( 1.0, epsSubst, 1.0, 1.0, incAngle(incWaveVec), k )*fe.rp(1.0, epsSub, 1.0, 1.0, angle, k )
    
def incAngle( incWaveVec ):
    k = np.sqrt( np.sum(incWaveVec**2) )
    return np.arccos( incWaveVec[2]/k )*180.0/np.pi

def qPar( incWaveVec ):
    return np.sqrt(incWaveVec[0]**2 + incWaveVec[1]**2)

def kScatteredPerp( incWaveVec, x, y, z ):
    k = np.sqrt( np.sum(incWaveVec**2) )
    return k*np.sqrt(x**2 + y**2)/np.sqrt(x**2+y**2+z**2)
    
def formFactorSphere( q_parallel, q_norm, Rsphere ):
    qR = np.sqrt( q_parallel**2 + q_norm**2 )*Rsphere
    qR[np.abs(qR) < 1E-5] = 1E-5
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
    k[1] = -kR*np.sin(grazingAngle*alpha_c)
    k[2] = kR*np.cos(grazingAngle*alpha_c)

    k_scattered = np.zeros(3)
    k_scattered[1] = kR*np.sin(grazingAngle*alpha_c)
    k_scattered[2] = kR*np.cos(grazingAngle*alpha_c)

    incAngle = 90.0 - grazingAngle*alpha_c*180.0/np.pi
    if ( polarisation == "TE" ):
        r = fe.rs( epsInc, epsSub, mu, mu, incAngle, kR )
    elif ( polarisation == "TM" ):
        r = fe.rp( epsInc, epsSub, mu, mu, incAngle, kR )
    else:
        print ("Unknown polarisation...")
        return
    
    z = 1000.0
    x = np.linspace( -z, z, 101 )
    y = np.linspace( -z, z, 101 )
    X, Y = np.meshgrid( x, y )
    f1 = Born1( k, X, Y, z, Rsphere )
    f2 = reflectedFromSubsrate( k_scattered, X, Y, z, Rsphere, epsSub )
    f3 = reflectedFromSubstrateAfter( k_scattered, X, Y, z, Rsphere, epsSub )
    f4 = reflectedFromSubsrateBeforeAndAfter( k_scattered, X, Y, z, Rsphere, epsSub )
    f = np.abs(f1 + f2 + f3 + f4)**2
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.contourf(X,Y,f,200, map="gist_heat")
    
    fname = "Figures/pattern.png"
    fig.savefig( fname, bbox_inches="tight" )
    print ("Figure written to %s"%(fname))

    f1 = Born1( k, x, 0.0, z, Rsphere )
    f2 = reflectedFromSubsrate( k_scattered, x, 0.0, z, Rsphere, epsSub )
    f3 = reflectedFromSubstrateAfter( k_scattered, x, 0.0, z, Rsphere, epsSub )
    f4 = reflectedFromSubsrateBeforeAndAfter( k_scattered, x, 0.0, z, Rsphere, epsSub )
    f = np.abs(f1 + f2 + f3 + f4)**2
     
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot( x, f1, 'k', label="Parallel" )
    ax.set_xlabel("Distance from center of screen")
    ax.legend(frameon=False)
    fname = "Figures/pattern1D.pdf"
    fig.savefig( fname, bbox_inches="tight" )
    print ("Figure written to %s"%(fname))

def main():
    grazingIncidence( 10.0, "TE", (1-1E-5+1j*1E-6)**2 )

if __name__ == "__main__":
    main() 
