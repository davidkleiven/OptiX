# OptiX
FDTD simulations of optical phenomena

# Dependencies
* HDF5 *-lhdf5* <https://www.hdfgroup.org/HDF5/doc/cpplus_RM/>
* Zlib *-lz* <http://zlib.net/>
* GNU Scientific Library (GSL) *-lgsl* <https://www.gnu.org/software/gsl/>
* Harminv *-lharminv* <http://ab-initio.mit.edu/wiki/index.php/Harminv>
* Lapack *-llapack* <http://www.netlib.org/lapack/>
* Cblas *-lcblas* <http://www.netlib.org/blas/>
* Atlas *-latlas* <http://math-atlas.sourceforge.net/>
* FFTW3 *-lfftw3* <http://www.fftw.org/>
* Jsoncpp *-ljsoncpp* <https://github.com/open-source-parsers/jsoncpp>
* MEEP *-lmeep* <http://ab-initio.mit.edu/wiki/index.php/Meep>  
* C math library *-lm*

# Python visualisation
The plots are generated using [matplotlib](http://matplotlib.org/) and [Mayavi](http://docs.enthought.com/mayavi/mayavi/).
When using Mayavi on Unix system errors concerning the use of wx occure (at least it did do me).
This was fixed by changing the backend to qt4 by typing

```bash
export ETS_TOOLKIT=qt4
```

in the command line. This may be added to the *.profile* for permanent change.

# Create Movie From Image Files
To create a movie from image files using *ffmpeg* type

```bash
ffmpeg -i visualize%d.png outfile.mp4
```

This assumes that the image files are named *visualize0.png*, *visualize1.png*, ..., *visualizeN.png*.

# SCUFF-EM and BEM
Notes on the running the boundary element code using the *SCUFF-EM* library can be found [here](FresnelBEM/Notes.md).
