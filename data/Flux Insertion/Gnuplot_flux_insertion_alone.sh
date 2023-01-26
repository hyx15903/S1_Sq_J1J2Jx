#set term postscript color
set term epslatex color colortext size 10cm,5cm standalone
set output "flux_insertion.tex"

a = 11
b = 2
d = 9
e = 4
f = 5
g = 6
h = 7
i = 8
c = 1.4
border_w = 4
z = 3.5

set multiplot layout 2,3 rowsfirst
#1
set ylabel "$ -ln(\\lambda _{i})$" offset 1
set yr [0:5]
set xr [-0.05:2.05]
set xlabel "$ \\theta /2\\pi $" offset 0,0.9
set lmargin at screen 0.1
set rmargin at screen 0.48
set tmargin at screen 0.92
set bmargin at screen 0.16
set border lw border_w
set xtics 0,0.5 offset 0,0.5
set ytics offset 0.5
set key samplen 0.5
set key spacing 0.5
set key at 2,0.9
set label "$ (a)$" at 0,5.3
set label "\\tiny $ L_{y}=8 $" at 0,0.4
plot "Flux_insertion.txt" index 0 u 1:(-log($2)) w lp title "\\tiny $ S_{z}=-2 $" pt a lc b ps c lw z,\
"Flux_insertion.txt" index 3 u 1:(-log($2)) w lp notitle pt a lc b ps c lw z,\
"Flux_insertion.txt" index 6 u 1:(-log($2)) w lp notitle pt a lc b ps c lw z,\
"Flux_insertion.txt" index 9 u 1:(-log($2)) w lp notitle pt a lc b ps c lw z,\
"Flux_insertion.txt" index 12 u 1:(-log($2)) w lp notitle pt a lc b ps c lw z,\
"Flux_insertion.txt" index 15 u 1:(-log($2)) w lp notitle pt a lc b ps c lw z,\
"Flux_insertion.txt" index 1 u 1:(-log($2)) w lp title "\\tiny $ S_{z}=-1 $" pt d lc e ps c lw z,\
"Flux_insertion.txt" index 4 u 1:(-log($2)) w lp notitle pt d lc e ps c lw z,\
"Flux_insertion.txt" index 7 u 1:(-log($2)) w lp notitle pt d lc e ps c lw z,\
"Flux_insertion.txt" index 10 u 1:(-log($2)) w lp notitle pt d lc e ps c lw z,\
"Flux_insertion.txt" index 13 u 1:(-log($2)) w lp notitle pt d lc e ps c lw z,\
"Flux_insertion.txt" index 16 u 1:(-log($2)) w lp notitle pt d lc e ps c lw z,\
"Flux_insertion.txt" index 2 u 1:(-log($2)) w lp title "\\tiny $ S_{z}=0 $" pt f lc g ps c lw z,\
"Flux_insertion.txt" index 5 u 1:(-log($2)) w lp notitle pt f lc g ps c lw z,\
"Flux_insertion.txt" index 8 u 1:(-log($2)) w lp notitle pt f lc g ps c lw z,\
"Flux_insertion.txt" index 11 u 1:(-log($2)) w lp notitle pt f lc g ps c lw z,\
"Flux_insertion.txt" index 14 u 1:(-log($2)) w lp notitle pt f lc g ps c lw z,\
"Flux_insertion.txt" index 17 u 1:(-log($2)) w lp notitle pt f lc g ps c lw z



#2
reset
set ylabel "$ -ln(\\lambda _{i})$" offset 1
set yr [0:5]
set xr [-0.05:1.05]
set xlabel "$ \\theta /2\\pi $" offset 0,0.9
set lmargin at screen 0.6
set rmargin at screen 0.98
set tmargin at screen 0.92
set bmargin at screen 0.16
set border lw border_w
set xtics 0,0.2 offset 0,0.5
set ytics offset 0.5
set key samplen 0.5
set key spacing 0.5
set key at 1,0.9
set label "$ (b)$" at 0,5.3
set label "\\tiny $ L_{y}=7 $" at 0,0.4
plot "Flux_insertion_Nx_2.txt" index 0 u 1:(-log($2)) w lp title "\\tiny $ S_{z}=-1 $" pt d lc e ps c,\
"Flux_insertion_Nx_2.txt" index 3 u 1:(-log($2)) w lp notitle pt d lc e ps c lw z,\
"Flux_insertion_Nx_2.txt" index 6 u 1:(-log($2)) w lp notitle pt d lc e ps c lw z,\
"Flux_insertion_Nx_2.txt" index 9 u 1:(-log($2)) w lp notitle pt d lc e ps c lw z,\
"Flux_insertion_Nx_2.txt" index 12 u 1:(-log($2)) w lp notitle pt d lc e ps c lw z,\
"Flux_insertion_Nx_2.txt" index 15 u 1:(-log($2)) w lp notitle pt d lc e ps c lw z,\
"Flux_insertion_Nx_2.txt" index 1 u 1:(-log($2)) w lp title "\\tiny $ S_{z}=0 $" pt f lc g ps c lw z,\
"Flux_insertion_Nx_2.txt" index 4 u 1:(-log($2)) w lp notitle pt f lc g ps c lw z,\
"Flux_insertion_Nx_2.txt" index 7 u 1:(-log($2)) w lp notitle pt f lc g ps c lw z,\
"Flux_insertion_Nx_2.txt" index 10 u 1:(-log($2)) w lp notitle pt f lc g ps c lw z,\
"Flux_insertion_Nx_2.txt" index 13 u 1:(-log($2)) w lp notitle pt f lc g ps c lw z,\
"Flux_insertion_Nx_2.txt" index 16 u 1:(-log($2)) w lp notitle pt f lc g ps c lw z,\
"Flux_insertion_Nx_2.txt" index 2 u 1:(-log($2)) w lp title "\\tiny $ S_{z}=1 $" pt h lc rgb "brown" ps c lw z,\
"Flux_insertion_Nx_2.txt" index 5 u 1:(-log($2)) w lp notitle pt h lc rgb "brown" ps c lw z,\
"Flux_insertion_Nx_2.txt" index 8 u 1:(-log($2)) w lp notitle pt h lc rgb "brown" ps c lw z,\
"Flux_insertion_Nx_2.txt" index 11 u 1:(-log($2)) w lp notitle pt h lc rgb "brown" ps c lw z,\
"Flux_insertion_Nx_2.txt" index 14 u 1:(-log($2)) w lp notitle pt h lc rgb "brown" ps c lw z,\
"Flux_insertion_Nx_2.txt" index 17 u 1:(-log($2)) w lp notitle pt h lc rgb "brown" ps c lw z




unset multiplot
