import numpy as np
def rs( eps1, eps2, mu1, mu2, angle, k ):
    n1 = np.sqrt(eps1*mu1)
    n2 = np.sqrt(eps2*mu2)
    kxi = k*np.sin(angle*np.pi/180.0)
    kzi = k*np.cos(angle*np.pi/180.0)
    kzt = np.sqrt( k**2 - (n1*kxi/n2)**2 +0j )
    a1 = np.sqrt( eps1/mu1 )
    a2 = np.sqrt( eps2/mu2 )
    print kzt, k
    return ( a1*kzi - a2*kzt )/( a1*kzi + a2*kzt )

def rp( eps1, eps2, mu1, mu2, angle, k ):
    a1 = np.sqrt( eps1/mu1 )
    a2 = np.sqrt( eps2/mu2 )
    kxi = k*np.sin(angle*np.pi/180.0)
    kzi = k*np.cos(angle*np.pi/180.0)
    n1 = np.sqrt(eps1*mu1)
    n2 = np.sqrt(eps2*mu2)
    kzt = np.sqrt( k**2 - (n1*kxi/n2)**2 +0j )
    return ( a2*kzi - a1*kzt )/( a1*kzt + a2*kzi )

def ts( eps1, eps2, mu1, mu2, angle, k ):
    a1 = np.sqrt( eps1/mu1 )
    a2 = np.sqrt( eps2/mu2 )
    kxi = k*np.sin(angle*np.pi/180.0)
    kzi = k*np.cos(angle*np.pi/180.0)
    n1 = np.sqrt(eps1*mu1)
    n2 = np.sqrt(eps2*mu2)
    kzt = np.sqrt( k**2 - (n1*kxi/n2)**2 )
    return 2.0*kzi*a1/( a1*kzi + a2*kzt )

 
def tp( eps1, eps2, mu1, mu2, angle, k ):
    a1 = np.sqrt( eps1/mu1 )
    a2 = np.sqrt( eps2/mu2 )
    kxi = k*np.sin(angle*np.pi/180.0)
    kzi = k*np.cos(angle*np.pi/180.0)
    n1 = np.sqrt(eps1*mu1)
    n2 = np.sqrt(eps2*mu2)
    kzt = np.sqrt( k**2 - (a1*kxi/a2)**2 )
    return 2.0*kzi*a1/( a1*kzt + a2*kzi )
    

# Reflectance and transmittance
def Rs( eps1, eps2, mu1, mu2, kzi, k ):
    #print rs(eps1,eps2,mu1,mu2,kzi,k)
    return np.abs( rs(eps1,eps2,mu1,mu2,kzi,k) )**2

def Rp( eps1, eps2, mu1, mu2, kzi, k ):
    return np.abs( rp(eps1,eps2,mu1,mu2,kzi,k) )**2

def Ts( eps1, eps2, mu1, mu2, kzi, k ):
    return np.abs( rp(eps1,eps2,mu1,mu2,kzi,k) )**2
