set terminal pdfcairo
set output "Figures/planeReflection/ampRatio.pdf"

FNAME = "dataPlane/transmittance.csv"
set datafile separator ","

n1 = 1.0
n2 = 1.5
PI=acos(-1.0)

cost(x) = sqrt(1.0 - (n1*sin(x*PI/180.0)/n2)**2)
# Exact values
ts(x) = 2.0*n1*cos(x*PI/180.0)/( n1*cos(x*PI/180.0) + n2*cost(x) )

set xlabel "Incident angle"
set ylabel "Amplitude ratio"
set xrange[0:90]
set yrange[0:1]
plot FNAME using 1:(sqrt($2)) with points pt 6 lc "black" notitle,\
ts(x) with lines lc "black" notitle
