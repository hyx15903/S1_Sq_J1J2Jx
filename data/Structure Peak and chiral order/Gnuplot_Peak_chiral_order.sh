set term epslatex color size 12cm,12cm standalone
#set term postscript color
set output "Peak_chiral_order.tex"

xrange = 0.45
s=1.3
border_w = 4
line_w = 3

set multiplot layout 2,3 rowsfirst
#1
set ylabel "$m^{2}(0,\\pi)$" offset 2.5
set xlabel "$ 1/L_{y}$" offset 0,0.7
set lmargin at screen 0.1
set rmargin at screen 0.48
set tmargin at screen 0.96
set bmargin at screen 0.58
set border lw border_w
set key samplen 0.5
set key spacing 0.7
set key at 0.17,0.386
set xr [0:0.25]
set yr [0:0.4]
set xtics 0,0.05 offset 0,0.5
set ytics 0,0.05 offset 0.7
set format x '%g'
set format y '%g'

set label "(a)" at 0.005,0.42

f0(x) = a0*x*x+b0*x*x*x*x+c0
fit f0(x) "Peak_value_j2_0.5_Nx_2.txt" index 0 using (1/$1):($2/$1/$1) via a0,b0,c0
f1(x) = a1*x*x+b1*x*x*x*x+c1
fit f1(x) "Peak_value_j2_0.5_Nx_2.txt" index 1 using (1/$1):($2/$1/$1) via a1,b1,c1
f2(x) = a2*x*x+b2*x*x*x*x+c2
fit f2(x) "Peak_value_j2_0.5_Nx_2.txt" index 2 using (1/$1):($2/$1/$1) via a2,b2,c2
f3(x) = a3*x*x+b3*x*x*x*x+c3
fit f3(x) "Peak_value_j2_0.5_Nx_2.txt" index 3 using (1/$1):($2/$1/$1) via a3,b3,c3
g2(x) = d2*x*x*x*x+e2*x*x+f2
fit g2(x) "scaling_order_j2_0.5_jchi_0.3.txt" index 0 u (1/$1):2 via d2,e2,f2


plot "Peak_value_j2_0.5_Nx_2.txt" index 0 using (1/$1):($2/$1/$1) w p title "\\tiny $ J_{2}=0.5, J_{\\chi }=0.19$" lw line_w lc rgb "black" pt 5 ps s,\
f0(x) notitle lc rgb "black" lw line_w,\
"Peak_value_j2_0.5_Nx_2.txt" index 1 using (1/$1):($2/$1/$1) w p title "\\tiny $ J_{2}=0.5, J_{\\chi }=0.2$" lw line_w lc rgb "blue" pt 7 ps s,\
f1(x) notitle lc rgb "blue" lw line_w,\
"Peak_value_j2_0.5_Nx_2.txt" index 2 using (1/$1):($2/$1/$1) w p title "\\tiny $ J_{2}=0.5, J_{\\chi }=0.21$" lw line_w lc rgb "red" pt 9 ps s,\
f2(x) notitle lc rgb "red" lw line_w,\
"Peak_value_j2_0.5_Nx_2.txt" index 3 using (1/$1):($2/$1/$1) w p title "\\tiny $ J_{2}=0.5, J_{\\chi }=0.22$" lw line_w lc rgb "dark-goldenrod" pt 13 ps s,\
f3(x) notitle lc rgb "dark-goldenrod" lw line_w,\
"scaling_order_j2_0.5_jchi_0.3.txt" index 0 u (1/$1):2 w p title "\\tiny $ J_{2}=0.5,J_{\\chi }=0.3 $" lw line_w pt 11 lc 1 ps s,\
g2(x) notitle lc 1 lw line_w


#2
reset
set ylabel "$m^{2}(\\pi ,\\pi )$" offset 2.4
set xlabel "$ 1/L_{y}$" offset 0,0.7
set lmargin at screen 0.59
set rmargin at screen 0.97
set tmargin at screen 0.96
set bmargin at screen 0.58
set border lw border_w
set key samplen 0.5
set key spacing 0.7
set key at 0.17,0.386
set xr [0:0.25]
set yr [0:0.4]
set xtics 0,0.05 offset 0,0.5
set ytics 0,0.05 offset 0.7
set format x '%g'
set format y '%g'
set label "(b)" at 0.005,0.42

g0(x) = d0*x*x+e0*x+f0
fit g0(x) "Peak_value_j2_0.5_Nx_2.txt" index 0 using (1/$1):($4/$1/$1) via d0,e0,f0
g1(x) = d1*x*x+e1*x+f1
fit g1(x) "Peak_value_j2_0.5_Nx_2.txt" index 1 using (1/$1):($4/$1/$1) via d1,e1,f1
g2(x) = d2*x*x+e2*x+f2
fit g2(x) "Peak_value_j2_0.5_Nx_2.txt" index 2 using (1/$1):($4/$1/$1) via d2,e2,f2
g3(x) = d3*x*x+e3*x+f3
fit g3(x) "Peak_value_j2_0.5_Nx_2.txt" index 3 using (1/$1):($4/$1/$1) via d3,e3,f3


plot "Peak_value_j2_0.5_Nx_2.txt" index 0 using (1/$1):($4/$1/$1) w p title "\\tiny $ J_{2}=0.5, J_{\\chi }=0.19$" lw line_w lc rgb "black" pt 4 ps s,\
g0(x) notitle lc rgb "black" lw line_w,\
"Peak_value_j2_0.5_Nx_2.txt" index 1 using (1/$1):($4/$1/$1) w p title "\\tiny $ J_{2}=0.5, J_{\\chi }=0.2$" lw line_w lc rgb "blue" pt 6 ps s,\
g1(x) notitle lc rgb "blue" lw line_w,\
"Peak_value_j2_0.5_Nx_2.txt" index 2 using (1/$1):($4/$1/$1) w p title "\\tiny $ J_{2}=0.5, J_{\\chi }=0.21$" lw line_w lc rgb "red" pt 8 ps s,\
g2(x) notitle lc rgb "red" lw line_w,\
"Peak_value_j2_0.5_Nx_2.txt" index 3 using (1/$1):($4/$1/$1) w p title "\\tiny $ J_{2}=0.5, J_{\\chi }=0.22$" lw line_w lc rgb "dark-goldenrod" pt 12 ps s,\
g3(x) notitle lc rgb "dark-goldenrod" lw line_w


#3
reset
set ylabel "$ <\\chi >$" offset 2
set xlabel "$ J_{\\chi }$" offset 0,0.7
set lmargin at screen 0.1
set rmargin at screen 0.48
set tmargin at screen 0.46
set bmargin at screen 0.08
set border lw border_w
set key samplen 0.5
set key spacing 0.5
set key right bottom
set yr [0:0.6]
set xr [0:xrange]
set xtics 0,0.1 offset 0,0.5
set ytics 0,0.1 offset 0.7
set format x '%g'
set format y '%g'
set label "(c)" at 0.006,0.63
set arrow 1 from 0.21,0 to 0.21,0.25 nohead dt 3 lc -1 lw 4.5
set arrow 3 from 0.21,0.524 to 0.21,0.6 nohead dt 3 lc 0 lw 4.5
set arrow 2 from 0.39,0 to 0.39,0.6 nohead dt 5 lc rgb "black" lw 4

plot "J20.5.txt" index 1 u 1:(abs($2)/4) w lp title "\\tiny $ J_{2}=0.5$" lw line_w pt 2 lc 2 ps s

#3 inset
reset
unset ylabel
set xlabel "\\tiny $ J_{\\chi }$" offset 0,1.5
set lmargin at screen 0.14
set rmargin at screen 0.29
set tmargin at screen 0.41
set bmargin at screen 0.26
set border lw border_w
#set yr [0:0.7]
set key samplen 0.5
set key spacing 0.5
set xr [0:xrange]
set xtics 0,0.1 offset 0,0.7
set ytics 0,2 offset 0.7
set format x '\tiny %g'
set format y '\tiny %g'
x0=NaN
y0=NaN

plot "J20.5.txt" index 1 u (dx=$1-x0,x0=$1,$1-dx/2):(dy=(abs($2)/4)-y0,y0=(abs($2)/4),dy/dx) w lp title "\\tiny $ \\frac{d<\\chi >}{dJ_{\\chi }} $" lw 2 pt 2 lc 2 ps s

#4
reset
set ylabel "$ <\\chi >$" offset 2
set xlabel "$ J_{2}$" offset 0,0.7
set lmargin at screen 0.59
set rmargin at screen 0.97
set tmargin at screen 0.46
set bmargin at screen 0.08
set border lw border_w
set key samplen 0.5
set key spacing 0.5
set key at 0.485,0.33
set yr [0:0.35]
set xr [0.4:0.6]
set xtics 0.4,0.05 offset 0,0.5
set ytics 0,0.1 offset 0.7
set format x '%g'
set format y '%g'
set label "(d)" at 0.405,0.365
set arrow 4 from 0.475,0 to 0.475,0.35 nohead dt 3 lc -1 lw 4.5
set arrow 5 from 0.545,0 to 0.545,0.04 nohead dt 5 lc 0 lw 4.5
set arrow 6 from 0.545,0.19 to 0.545,0.35 nohead dt 5 lc rgb "black" lw 4


plot "chiral_order_J2.txt" index 0 u 1:(abs($2)/4) w lp title "\\tiny $ J_{\\chi}=0.25$" lw line_w pt 1 lc 2 ps s

#4 inset
reset
unset ylabel
set xlabel "\\tiny $ J_{\\chi }$" offset 0.3,1.7
set lmargin at screen 0.8
set rmargin at screen 0.95
set tmargin at screen 0.29
set bmargin at screen 0.14
set border lw border_w
#set yr [0:0.7]
set key samplen 0.5
set key spacing 0.5
set xr [0.4:0.6]
set yr [-1:7]
set xtics 0.4,0.1 offset 0,0.7
set ytics 0,2 offset 0.7
set format x '\tiny %g'
set format y '\tiny %g'
#set arrow 7 from 0.4,0 to 0.6,0 nohead dt 6 lc 0 lw 0.7
x0=NaN
y0=NaN

plot "chiral_order_J2.txt" index 0 u (dx=$1-x0,x0=$1,$1-dx/2):(dy=(abs($2)/4)-y0,y0=(abs($2)/4),dy/dx) w lp title "\\tiny $ \\frac{d<\\chi >}{dJ_{2 }} $" lw 2 pt 1 lc 2 ps s


unset multiplot
