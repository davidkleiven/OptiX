import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 28
mpl.rcParams["axes.unicode_minus"] = False

delta = 4.14E-5
widths = [50.0, 100.0, 200.0, 500.0, 1000.0] # nano meter
colors = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99"]
k = 40.05 # in nm^{-1}

def main():
    assert( len(widths) == len(colors) )
    deltaKK = delta*k*k # This is a constant independent of k in the high frequency approximation (delta ~ k^{-2})

    print ( "Drude proportionality: %.2E"%(deltaKK))

    fig = plt.figure()
    axSym = fig.add_subplot(1,2,1)
    fig.suptitle( "\$ \frac{Ne^2}{mc^2} = %.2f nm^{-2}\$"%(2.0*deltaKK) )

    axAsym = fig.add_subplot(1,2,2)
    xmax = k*widths[-1]/100.0
    start = 0.0
    counter = 0
    while ( start**2 < xmax ):
        start = counter*np.pi/2.0+1E-12
        end = (counter+1)*np.pi/2.0 - 1E-12
        xtemp = np.linspace(start**2, end**2, 101)
        symLhs = np.sqrt(xtemp)*np.tan(np.sqrt(xtemp) )
        aSymLhs = np.sqrt(xtemp)/np.tan(np.sqrt(xtemp))
        axSym.plot( xtemp, symLhs, color="black")
        axAsym.plot( xtemp, aSymLhs, color="black")
        start += ((np.pi/2.0)**2+1E-12)
        counter += 1

    x = np.linspace(0.0, 1.3*xmax, 1001)
    for i in range( 0, len(widths) ):
        qwHalf = np.sqrt( (k*widths[i]/2.0)**2 - x**2)

        rhs = 0.5*widths[i]*np.sqrt( 2.0*deltaKK - (2.0/widths[i]**2 *x) )
        axSym.plot( x,  rhs, color=colors[i], label=widths[i] )
        axAsym.plot( x, rhs, color=colors[i], label=widths[i] )
    ymax = 0.8*widths[-1]*np.sqrt(2.0*deltaKK)
    axSym.set_ylim(0.0, ymax)
    axSym.set_xlabel( "\$\frac{qw}{2}\$" )
    axAsym.set_xlabel( "\$\frac{qw}{2}\$" )
    axAsym.set_ylim(0.0, ymax)
    plt.show()

if __name__ == "__main__":
    main()
