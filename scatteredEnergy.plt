set terminal epslatex size 3.5,2.62
set output "Figures/energyScattered.tex"

FNAME1="dataPlane/MultInc45/WithEps/transmittedFluxNorm.csv"
FNAME2="dataPlane/MultInc20/WithEps/transmittedFluxNorm.csv"
FNAME3="dataPlane/MultInc75/WithEps/transmittedFluxNorm.csv"
FNAME4="dataPlane/MultInc85/WithEps/transmittedFluxNorm.csv"
FNAME1p="dataPlane/MultInc45p/WithEps/transmittedFluxNorm.csv"
FNAME2p="dataPlane/MultInc20p/WithEps/transmittedFluxNorm.csv"
FNAME3p="dataPlane/MultInc75p/WithEps/transmittedFluxNorm.csv"
FNAME4p="dataPlane/MultInc85p/WithEps/transmittedFluxNorm.csv"
set datafile separator ","
n1 = 1.0
n2 = 1.5

PI=acos(-1.0)
cosTx(x) = sqrt(1.0 - (n1*sin(x*PI/180.0)/n2)**2)
Rs(x) = (( n1*cos(x*PI/180.0) - n2*cosTx(x) )/( n1*cos(x*PI/180.0) + n2*cosTx(x) ))**2
Rp(x) = (( n1*cosTx(x) - n2*cos(x*PI/180.0) )/( n2*cos(x*PI/180.0) + n1*cosTx(x) ))**2

Ts(x) = 1.0 - Rs(x)
Tp(x) = 1.0 - Rp(x)
set xlabel "Incident angle (deg)"
set ylabel "Scattered power"
set yrange[0:1]
set xrange[0:90]
set key center left
plot FNAME1 every 10 using 1:2 with points pt 6 ps 0.3 lc "black" title "$T_s$",\
Ts(x) with lines lt 1 lc "black" lw 4 notitle,\
FNAME2 every 10 using 1:2 with points pt 6 ps 0.3 lc "black" notitle,\
FNAME3 every 1 using 1:2 with points pt 6 ps 0.3 lc "black" notitle,\
FNAME4 every 1 using 1:2 with points pt 6 ps 0.3 lc "black" notitle,\
Tp(x) with lines lt 3 lc "black" lw 4 notitle,\
FNAME1p every 10 using 1:2 with points pt 8 ps 0.4 lc "black" title "$T_p$",\
FNAME2p every 10 using 1:2 with points pt 8 ps 0.4 lc "black" notitle,\
FNAME3p every 1 using 1:2 with points pt 8 ps 0.4 lc "black" notitle,\
FNAME4p every 1 using 1:2 with points pt 8 ps 0.4 lc "black" notitle,\
FNAME1 every 10 using 1:3 with points pt 10 ps 0.4 lc "black" title "$R_s$",\
FNAME2 every 10 using 1:3 with points pt 10 ps 0.4 lc "black" notitle,\
FNAME3 every 1 using 1:3 with points pt 10 ps 0.4 lc "black" notitle,\
FNAME4 every 1 using 1:3 with points pt 10 ps 0.4 lc "black" notitle,\
Rs(x) with lines lt 1 lc "black" lw 4 notitle,\
FNAME1p every 10 using 1:3 with points pt 4 ps 0.3 lc "black" title "$R_p$",\
FNAME2p every 10 using 1:3 with points pt 4 ps 0.3 lc "black" notitle,\
FNAME3p every 1 using 1:3 with points pt 4 ps 0.3 lc "black" notitle,\
FNAME4p every 1 using 1:3 with points pt 4 ps 0.3 lc "black" notitle,\
Rp(x) with lines lt 3 lc "black" lw 4 notitle

set terminal pdfcairo dashed
set output "Figures/energyScattered.pdf"
replot
