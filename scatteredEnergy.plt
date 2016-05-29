set terminal pdfcairo
set output "Figures/energyScattered.pdf"

FNAME1="dataPlane/MultInc45/WithEps/transmittedFluxNorm.csv"
FNAME2="dataPlane/MultInc20/WithEps/transmittedFluxNorm.csv"
FNAME3="dataPlane/MultInc85/WithEps/transmittedFluxNorm.csv"
set datafile separator ","
n1 = 1.0
n2 = 1.5

PI=acos(-1.0)
cosTx(x) = sqrt(1.0 - (n1*sin(x*PI/180.0)/n2)**2)
Rs(x) = (( n1*cos(x*PI/180.0) - n2*cosTx(x) )/( n1*cos(x*PI/180.0) + n2*cosTx(x) ))**2
Ts(x) = 1.0 - Rs(x)
set xlabel "Incident angle (deg)"
set ylabel "Scattered power"
set yrange[0:1]
set xrange[0:90]
plot FNAME1 every 10 using 1:2 with points pt 6 ps 0.3 lc "black" notitle,\
Ts(x) with lines lc "black" lw 4 title "s",\
FNAME2 every 10 using 1:2 with points pt 6 ps 0.3 lc "black" notitle,\
FNAME3 every 10 using 1:2 with points pt 6 ps 0.3 lc "black" notitle

