set terminal pngcairo
set output ofile

# ofile and ifile are command line arguments
set datafile separator ","

set xlabel "Time (au)"
set ylabel "Amplitude"
plot ifile using 1:2 with lines lc "black" notitle
