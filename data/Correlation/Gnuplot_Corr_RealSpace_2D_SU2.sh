set term epslatex color size 10cm,5cm standalone
#set term postscript color
set output "Correlation_2D_SU2.tex"

##set lmargin at screen 0.1
##set rmargin at screen 1.0
##set key outside right
set multiplot layout 2,1 rowsfirst

border_w = 4

#1
set xlabel "$ x $" offset 0,0.5
set ylabel "$ y $" offset 1.3
set lmargin at screen 0.07
set rmargin at screen 0.47
set tmargin at screen 0.9
set bmargin at screen 0.1
set border lw border_w
set xr [0.5:8.5]
set yr [0.5:8.5]
set xtics 1,1
set ytics 1,1
unset xtics
unset ytics
dx=1  #grid spacing in x
set for [i=1:8] arrow from i*dx,1 to i*dx,8 nohead back ls 2
dy=1  #grid spacing in y
set for [i=1:8] arrow from 1,i*dy to 8,i*dy nohead back ls 2

set object circle at 1,4 size scr 0.013 fs transparent fc rgb "black" lc rgb "black" front
set style textbox opaque noborder fillcolor rgb "white"
#set label "-0.08" at 1.8,3.8 front boxed

set label "(a)" at 0.5,9
#set label "0.0006" at 4.8,3.8 front boxed
set label "1" at 0.6,1
set label "2" at 0.6,2
set label "3" at 0.6,3

plot "spin_cor_J2_05_Jx_045_middle.dat" index 0 u ((int($2)-((int($2)-64-1)%8+1)-64)/8+1):((int($2)-64-1)%8+1):($3>0?(abs($3)):0) '          %lf          %lf (%lf,%lf)' with circle lc rgb "blue" fill solid border notitle,\
"spin_cor_J2_05_Jx_045_middle.dat" index 0 u ((int($2)-((int($2)-64-1)%8+1)-64)/8+1):((int($2)-64-1)%8+1):($3<0?(abs($3)):0) '          %lf          %lf (%lf,%lf)' with circle lc rgb "red" fs transparent solid 0.2 border notitle 

#2
reset
set lmargin at screen 0.56
set rmargin at screen 0.96
set tmargin at screen 0.9
set bmargin at screen 0.1
set border lw border_w
set xr [0.5:8.5]
set yr [0.5:8.5]
set xtics 1,1
set ytics 1,1
set xlabel "$ x $" offset 0,0.5
set ylabel "$ y $" offset 1.3
set label "1" at 0.6,1
set label "2" at 0.6,2
set label "3" at 0.6,3

unset xtics
unset ytics
dx=1  #grid spacing in x
set for [i=1:8] arrow from i*dx,1 to i*dx,8 nohead back ls 2
dy=1  #grid spacing in y
set for [i=1:8] arrow from 1,i*dy to 8,i*dy nohead back ls 2
set object circle at 1,5 size scr 0.013 fs transparent fc rgb "black" lc rgb "black" front
set style textbox opaque noborder fillcolor rgb "white"

set label "(b)" at 0.5,9
#set label "0.04" at 4.8,3.8 front boxed

plot "spin_cor_J2_05_Jx_045_middle.dat" index 1 u ((int($2)-((int($2)-64-1)%8+1)-64)/8+1):((int($2)-64-1)%8+1):($3>0?(abs($3)):0) '          %lf          %lf (%lf,%lf)' with circle lc rgb "blue" fs solid border notitle,\
"spin_cor_J2_05_Jx_045_middle.dat" index 1 u ((int($2)-((int($2)-64-1)%8+1)-64)/8+1):((int($2)-64-1)%8+1):($3<0?(abs($3)):0) '          %lf          %lf (%lf,%lf)' with circle lc rgb "red" fs transparent solid 0.2 border notitle 


unset multiplot
