reset

set xlabel 'cm-3' font ',10' 
set ylabel 'erg s-1' font ',10' 

plot "t.data" using 1:2 with lines linewidth 6 title "Density"

pause -1
