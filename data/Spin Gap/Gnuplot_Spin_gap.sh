#set term postscript color
set term epslatex color colortext size 13cm,6.5cm standalone
set output "Spin_gap.tex"

border_w = 4
line_w = 3

a = 11
b = 2
d = 9
e = 4
f = 5
g = 6
h = 7
i = 8
c = 1
set multiplot layout 1,2 rowsfirst
#1
set lmargin at screen 0.1
set rmargin at screen 0.48
set tmargin at screen 0.9
set bmargin at screen 0.14
set border lw border_w
set key samplen 0.6
set key spacing 0.8
set ylabel "$ \\Delta_{s}$ " offset 2
set xlabel "$ 1/L_{x} $" offset 0,0.9
set key right bottom
set xr [0:0.08]
set yr [0.3:0.9]
set label "(a)" at 0.004,0.93
set xtics 0,0.02 offset 0,0.5
set ytics offset 0.5
plot "Spin_gap_SM.txt" index 2 using (1/$1):($3-$2) w lp title "$ J_{2}=0.45, J_{\\chi }=0.3$" lw line_w lc 3 pt 7 ps 2,\
"Spin_gap_SM.txt" index 3 using (1/$1):($3-$2) w lp title "$ J_{2}=0.5, J_{\\chi }=0.3$" lw line_w lc 4 pt 5 ps 2



#2
reset
set lmargin at screen 0.59
set rmargin at screen 0.97
set tmargin at screen 0.9
set bmargin at screen 0.14
set border lw border_w
set key samplen 0.6
set key spacing 0.8
set ylabel "$ \\Delta_{s}$"  offset 2
set xlabel "$ 1/L_{y}$" offset 0,0.9
set label "(b)" at 0.01,1.15
set xtics 0,0.05 offset 0,0.5
set ytics 0,0.2 offset 0.5
#f(x) = a0 *x*x + b0*x + c0
#fit f(x) "Spin_gap_SM.txt" index 0 using (1/$1/$1/3):($3-$2) via a0,b0,c0
#(x) = a1 *x*x + b1*x + c1
#fit g(x) "Spin_gap_SM.txt" index 1 using (1/$1/$1/3):($3-$2) via a1,b1,c1
#g(x) = a0 *x*x + b0*x - 0.03
set xr [0:0.25]
set yr [0:1.1]
set key right bottom

plot "Spin_gap_SM.txt" index 0 using (1/$1):($3-$2) w lp title "$ J_{2}=0.45, J_{\\chi }=0.3$" lw line_w lc 3 pt 6 ps 2,\
"Spin_gap_SM.txt" index 1 using (1/$1):($3-$2) w lp title "$ J_{2}=0.5, J_{\\chi }=0.3$" lw line_w lc 4 pt 4 ps 2
#((x > 0.004) ? g(x) : 1/0) notitle lc 4 dt 2

unset multiplot







