set key at 4,0.4

set title font "Verdana,12" 
set xtics font "Verdana,10" 
set ytics font "Verdana,10" 

set xrange [0.0:0.8]
set xlabel 'x' font ',10' 
set ylabel 'Pressure' font ',10' 
plot "riemann.data" using 1:4 with lines linewidth 1
replot "test1.txt" using 1:4
pause -1
set xlabel 'x' font ',10' 
set ylabel 'Density' font ',10' 
plot "riemann.data" using 1:2 with lines linewidth 1
replot "test1.txt" using 1:2
pause -1
set xlabel 'x' font ',10' 
set ylabel 'Energy' font ',10' 
plot "riemann.data" using 1:5 with lines linewidth 1
replot "test1.txt" using 1:3
pause -1
set xlabel 'x' font ',10' 
set ylabel 'Velocity' font ',10' 
plot "riemann.data" using 1:3 with lines linewidth 1
replot "test1.txt" using 1:7
pause -1
