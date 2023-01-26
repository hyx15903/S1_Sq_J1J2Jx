set term epslatex color colortext size 14cm,7cm standalone
#set term postscript color
set output "M_hz_order.tex"

border_w = 4
line_w = 3.5

set multiplot layout 2,1 rowsfirst
#1
set ylabel "$m^{2}(0,\\pi)$" offset 2.7
set xlabel "$ \\frac{1}{M}$" offset 0,0.8
set yr [0.17:0.2]
set xr [0:0.00065]
set key right bottom
set key samplen 0.8
set key spacing 0.8

set label "(a)" at 0.00002,0.2014
set lmargin at screen 0.11
set rmargin at screen 0.49
set tmargin at screen 0.9
set bmargin at screen 0.14
set border lw border_w
set xtics 0,0.0002 offset 0,0.5
set ytics offset 0.7
#set ytics 0,1,8 offset 0.9
#set label "$ S_{z} = 0$" at 0.8,0.5
g0(x) = d0*x*x+e0*x+f0
fit g0(x) "ms_J2_05_Jx_03.dat" index 1 u (1/$1):($2/8/8) via d0,e0,f0
g1(x) = d1*x*x+e1*x+f1
fit g1(x) "ms_J2_05_Jx_03.dat" index 2 u (1/$1):($2/10/10) via d1,e1,f1

plot "ms_J2_05_Jx_03.dat" index 1 u (1/$1):($2/8/8) w p title "$ L_{y}=8$" pt 6 lc 6 ps 2,\
(x < 0.0004167? g0(x):1/0) notitle lc 6 lw line_w,\
"ms_J2_05_Jx_03.dat" index 2 u (1/$1):($2/10/10) w p title "$ L_{y}=10 $" pt 7 lc 7 ps 2,\
g1(x) notitle lc 7 lw line_w




#2
reset
k=0.5
a=8
c=0
s=1.7
set xlabel "$ x$" offset 0,1
set ylabel "$ <S^{z}_{x}>$" offset 2
#set title "$ J_{2}=0.5, $"
set label "$ (b)$" at 0.7,0.84

#set key at 19,1.2
set key samplen 0.8
set key spacing 0.8

set lmargin at screen 0.61
set rmargin at screen 0.99
set tmargin at screen 0.9
set bmargin at screen 0.14
set border lw border_w

set format x '%g'
set xtics offset 0,0.5
set ytics offset 0.7
set yr [0:0.8]
set xr [0:15]
#a(x) = a0 * x + b0
#fit a(x) "StatSCorrNx20Ny8j20.01jz0jk10omega0test.txt" u ($1/8):(log(abs($2))) every 16::16 via a0,b0

#b(x) = a1 * x + b1
#fit b(x) "StatSCorrNx20Ny8j20.01jz0jk10omega0test.txt" u (log($1/8)):(log(abs($2))) every 16::16 via a1,b1

#set yr [log(0.0001):log(1)]
#set ylabel "ln(<S_{x_{0}} S_{x_{0}+x}>)"

plot "Pining_hz.txt" index 0 u (($1-1)/a):(abs($2)) every a::c w lp title "$ h_{z}=0.5$" lw line_w pt 6 lc 6 ps s,\
"Pining_hz.txt" index 1 u (($1-1)/a):(abs($2)) every a::c w lp title "$ h_{z}=0.1$" lw line_w pt 8 lc 8 ps s,\
"Pining_hz.txt" index 2 u (($1-1)/a):(abs($2)) every a::c w lp title "$ h_{z}=0.05$" lw line_w pt 10 lc 10 ps s,\
"Pining_hz.txt" index 3 u (($1-1)/a):(abs($2)) every a::c w lp title "$ h_{z}=0.01$" lw line_w pt 12 lc 12 ps s

#"Pining_hz.txt" index 2 u ($1/a):(abs($2)) every a::c w lp title "$ h_{z}=0.05$" lw line_w pt 10 lc 10,\


unset multiplot



