set term epslatex color colortext size 12cm,12cm standalone
#set term postscript color
set output "Bond_E_VEE_Corr_ESK.tex"

border_w = 4
line_w = 4


set multiplot layout 2,2 rowsfirst
#1
##set lmargin at screen 0.1
##set rmargin at screen 1.0
##set key outside right
set lmargin at screen 0.12
set rmargin at screen 0.47
set tmargin at screen 0.95
set bmargin at screen 0.6
set border lw border_w
set key samplen 0.5
set key spacing 0.5
set format x '\tiny %g'
#set format y '\tiny %g'
set xtics 0,0.0001 offset 0,0.7
set ytics -1.65,0.01 offset 0.7
set xlabel " $\\frac{1}{M} $" offset 0,1
set ylabel " $ E_{0}/N $" offset 2.6
set xr [0:0.0004]
set yr [-1.65:-1.60]
set label " (a)" at 0.00002,-1.598

plot "bond_E0_VEE_CorrLength.txt" index 0 u (1/$1):2 w lp title "\\tiny infinite DMRG" lw line_w pt 5 lc 11 ps 2.2,\
"bond_E0_VEE_CorrLength.txt" index 1 u (1/$1):2 w lp title "\\tiny finite DMRG" lw line_w pt 6 lc 9 ps 1.3




#2
reset
set lmargin at screen 0.62
set rmargin at screen 0.97
set tmargin at screen 0.95
set bmargin at screen 0.6
set border lw border_w
set key samplen 0.5
set key spacing 0.5
set xtics 0,0.0001 offset 0,0.7
set ytics offset 0.7
set format x '\tiny %g'
#set format y '\tiny %g'
set xlabel "$\\frac{1}{M} $" offset 0,1
set ylabel "entanglement entropy " offset 2.3
set xr [0:0.0004]
set yr [2:5]
set label "(b)" at 0.00002,5.1

plot "bond_E0_VEE_CorrLength.txt" index 0 u (1/$1):3 w lp title "\\tiny infinite DMRG" lw line_w pt 5 lc 11 ps 2.2,\
"bond_E0_VEE_CorrLength.txt" index 1 u (1/$1):3 w lp title "\\tiny finite DMRG" lw line_w pt 6 lc 9 ps 1.3



#3
reset
set lmargin at screen 0.12
set rmargin at screen 0.47
set tmargin at screen 0.45
set bmargin at screen 0.1
set border lw border_w
set key right top
#set format x '%2.0t*10^{%L}'
#set format y '\tiny %g'
set label "(c)" at 0.00002,1.55

set xlabel "$\\frac{1}{M} $ " offset 0,1
set ylabel "correlation length" offset 2
set xr [0:0.0003]
set yr [0:1.5]
#set format x "$%s*10^{%T}$"
set xtics 0,0.0001 offset 0,0.5
set ytics offset 0.7

plot "bond_E0_VEE_CorrLength.txt" index 0 u (1/$1):($4/8) w lp notitle lw line_w pt 7 lc 8 ps 1.8


#4
reset
delta = 0.15
x0 = 4.768
set lmargin at screen 0.62
set rmargin at screen 0.97
set tmargin at screen 0.45
set bmargin at screen 0.1
set border lw border_w
#set format x '\tiny %g'
#set format y '\tiny %g'
set ylabel "$ -ln(\\lambda _{i})$" offset 1.3
set xlabel "$\\Delta k_{y}$" offset 0,1
set yr [0:8.4]
set xr [0-x0:6.5-x0]
set xtics -6,2 offset 0,0.5
set ytics offset 0.7
set key at 2.7-x0,9.9
set key samplen 0.5
set key spacing 0.5
set label "(d)" at 1.7-x0,8
plot "bond_ESK.txt" index 1 u (($3 > -2.736)? ($3 + 2.1054 + delta*3-x0) : ($3 + 9.0884 + delta*3-x0)):(-log($2)) w p title "\\tiny $ M = 6000$" pt 4 lc 7 ps 1.2,\
"bond_ESK.txt" index 2 u (($3 > -2.136)? ($3 + 2.8354 - delta*5-x0) : ($3 + 9.1184 - delta*5-x0)):(-log($2)) w p title "\\tiny $ M = 8000$" pt 6 lc 8 ps 1.2,\
"bond_ESK.txt" index 3 u (($3 > -0.02)? ($3+0.02+delta*2-x0) : ($3 + 6.303+delta*2-x0)):(-log($2)) w p title "\\tiny $ M = 10000$" pt 8 lc rgb "blue" ps 1.2


unset multiplot

