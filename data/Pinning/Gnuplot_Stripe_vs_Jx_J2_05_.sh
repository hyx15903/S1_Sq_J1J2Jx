set term epslatex color size 9cm,9cm standalone
#set term postscript color
set output "Stripe_order_Jchi.tex"

set lmargin at screen 0.18
set rmargin at screen 0.98
set tmargin at screen 0.9
set bmargin at screen 0.1


k=0.5
a=8
c=0
set xlabel "$ J_{\\chi }$" offset 0,1
set ylabel "stripe order" offset 1.5
set title "$ J_{2}=0.5, N_{y}=10$ iDMRG"
#set label "$ (b) J_{2}=0.5$" at 0.3,1.5
#set key at 19,1.2
set key samplen 0.8
set key spacing 0.8

set xtics offset 0,0.5
set ytics offset 0.7
#set yr [-1:1]
#set xr [1:20]

plot "Stripe_Jx_iDMRG.txt" index 0 u 1:2 w lp notitle lw 1 pt 7 lc 7 ps 1.4
