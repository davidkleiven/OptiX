set terminal pngcairo
set output ofile

# ofile and ifile are command line arguments
set datafile separator ","

set xlabel "Time (au)"
set ylabel "Amplitude"
plot ifile using 1:2 with lines lc "black" title "Real",\
'' using 1:3 with lines lt 3 lc "blue" title "Imag"
