#set term postscript color
set term epslatex color colortext size 20cm,7cm standalone 
set output "state_entan_spec.tex"

a = 2
b = 6
c = 1.2
d = 0.6
e = 0.2
x = 4.489
x4 = 4.489
kx = 0.2
ky = 1.26
border_w = 4
line_w = 2
#1
set multiplot layout 1,10 rowsfirst
set ylabel "$ -ln(\\lambda _{i})$" offset 1.7
set yr [0:10.4]
set xr [-0.5-x:6.5-x]
set xlabel "$\\Delta k_{y}$" offset kx,ky
set label "(a) vacuum sector" at -0.5-x,10.65
set lmargin at screen 0.05
set rmargin at screen 0.13
set tmargin at screen 0.95
set bmargin at screen 0.1
set border lw border_w
set xtics -4,2,0 offset 0,0.5
set ytics offset 0.8
set label "\\tiny $ S_{z} = -1$" at 0.4-x,0.4
set label "1" at 3.803 - d-x,2.5 tc rgb "red"
set label "2" at 2.818 - d-x,2.9-e tc rgb "red"
set arrow from 2.918-x,-log(0.03803599553) to 2.918-x,-log(0.007708199) nohead lc rgb 'red' lw line_w
set label "4" at 2.1326 - d-x,3.5-e tc rgb 'red'
set arrow from 2.1326-x,-log(0.021772568146) to 2.1326-x,-log(0.0028469324808) nohead lc rgb 'red' lw line_w
set label "7" at 1.3473 - d-x,5.1-e tc rgb 'red'
set arrow from 1.3473-x,-log(0.004129401) to 1.3473-x,-log(0.00119626484798) nohead lc rgb 'red' lw line_w
set label "13" at 0.5619 - 2*d-x,5.7-e tc rgb 'red'
set arrow from 0.5619-x,-log(0.00244449342381402) to 0.5619-x,-log(3.7868382455748424e-05) nohead lc rgb 'red' lw line_w
plot "spectrum.txt" index 0 u (($3 > -0.02)? ($3 -x4) : ($3 + 6.283 -x4)):(-log($2)) w p notitle pt a lc b ps c

unset arrow
unset label
unset ylabel
set lmargin at screen 0.13
set rmargin at screen 0.21
set xlabel "$\\Delta  k_{y}$" offset kx,ky
unset ytics
set xtics -4,2,0 offset 0,0.5
set label "\\tiny $ S_{z} = 0$" at 1-x,0.4
set label "1" at 4.3887 - d-x,0.7 tc rgb "red"
set label "1" at 3.603 - d-x,2.5-e tc rgb "red"
set label "3" at 2.818 - d-x,2.9-e tc rgb "red"
set arrow from 2.918-x,-log(0.03803599553) to 2.918-x,-log(0.007708199) nohead lc rgb 'red' lw line_w
set label "5" at 2.1326 - d-x,3.5-e tc rgb 'red'
set arrow from 2.1326-x,-log(0.021772568146) to 2.1326-x,-log(0.0028469324808) nohead lc rgb 'red' lw line_w
set label "10" at 1.3473 - 2*d-x,4.8-e tc rgb 'red'
set arrow from 1.3473-x,-log(0.0060287329955) to 1.3473-x,-log(0.00098071501693) nohead lc rgb 'red' lw line_w
set label "16" at 0.5619 - 2*d-x,5.7-e tc rgb 'red'
set arrow from 0.5619-x,-log(0.00257791419281932) to 0.5619-x,-log(3.7868382455748424e-05) nohead lc rgb 'red' lw line_w
plot "spectrum.txt" index 1 u (($3 > -0.02)? ($3-x4) : ($3 + 6.283-x4)):(-log($2)) w p notitle pt a lc b ps c

unset arrow
unset label
set lmargin at screen 0.21
set rmargin at screen 0.29
set xlabel "$\\Delta  k_{y}$" offset kx,ky
unset ytics
set label "\\tiny $ S_{z} = 1$" at 1-x,0.4
set label "1" at 3.603 - d-x,2.5 tc rgb "red"
set label "2" at 2.818 - d-x,2.9-e tc rgb "red"
set arrow from 2.918-x,-log(0.03803599553) to 2.918-x,-log(0.007708199) nohead lc rgb 'red' lw line_w
set label "4" at 2.1326 - d-x,3.5-e tc rgb 'red'
set arrow from 2.1326-x,-log(0.021772568146) to 2.1326-x,-log(0.0028469324808) nohead lc rgb 'red' lw line_w
set label "7" at 1.3473 - d-x,5.1-e tc rgb 'red'
set arrow from 1.3473-x,-log(0.004129401) to 1.3473-x,-log(0.00119626484798) nohead lc rgb 'red' lw line_w
set label "13" at 0.5619 - 2*d-x,5.7-e tc rgb 'red'
set arrow from 0.5619-x,-log(0.00244449342381402) to 0.5619-x,-log(3.7868382455748424e-05) nohead lc rgb 'red' lw line_w
plot "spectrum.txt" index 2 u (($3 > -0.02)? ($3-x4) : ($3 + 6.283-x4)):(-log($2)) w p notitle pt a lc b ps c




x1 = 4.088
#2
reset
set lmargin at screen 0.36
set rmargin at screen 0.44
set tmargin at screen 0.95
set bmargin at screen 0.1
set border lw border_w
set xtics -4,2,0 offset 0,0.5
set ytics offset 0.8
set ylabel "$ -ln(\\lambda _{i})$" offset 1.7
set yr [0:10.4]
set xr [-0.5-x1:6.5-x1]
set xlabel "$\\Delta  k_{y}$" offset kx,ky
set label "(b) Ising anyon sector" at -0.5-x1,10.65
set label "\\tiny $ S_{z} = -1$" at 0.4-x1,0.4
set label "1" at 3.19-x1 - d,3 tc rgb "red"
set label "2" at 2.293 - d-x1,3.7-e tc rgb "red"
set arrow from 2.293-x1,-log(0.0169561123386) to 2.293-x1,-log(0.01019284373447120) nohead lc rgb 'red' lw line_w
set label "4" at 1.396 - d-x1,4.4-e tc rgb 'red'
set arrow from 1.396-x1,-log(0.008411066855) to 1.396-x1,-log(0.0023518969601) nohead lc rgb 'red' lw line_w
set label "8" at 0.498 - d-x1,5.5-e tc rgb 'red'
set arrow from 0.498-x1,-log(0.00294471015986) to 0.498-x1,-log(0.000364332075696) nohead lc rgb 'red' lw line_w
plot "spectrum.txt" index 3 u (($3 > 2.84)? ($3 - 3.23-x1) : ($3 + 3.053-x1)):(-log($2)) w p notitle pt a lc b ps c

unset arrow
unset label
unset ylabel
set lmargin at screen 0.44
set rmargin at screen 0.52
set xlabel "$\\Delta  k_{y}$" offset kx,ky
unset ytics
set label "\\tiny $ S_{z} = 0$" at 1-x1,0.4
set xtics -4,2,0 offset 0,0.5
set label "1" at 4.0878 - d-x1,1.5 tc rgb "red"
set label "2" at 3.19 - d-x1,2.8-e tc rgb "red"
set arrow from 3.19-x1,-log(0.04207561470370) to 3.19-x1,-log(0.033058813893777) nohead lc rgb 'red' lw line_w
set label "4" at 2.293 - d-x1,3.3-e tc rgb "red"
set arrow from 2.293-x1,-log(0.0234398254893215) to 2.293-x1,-log(0.0101927995217407) nohead lc rgb 'red' lw line_w
set label "8" at 1.396 - d-x1,4.4-e tc rgb 'red'
set arrow from 1.396-x1,-log(0.008411066855) to 1.396-x1,-log(0.0023518969601) nohead lc rgb 'red' lw line_w
set label "14" at 0.498 - 2*d-x1,5.5-e tc rgb 'red'
set arrow from 0.498-x1,-log(0.002944654050485286) to 0.498-x1,-log(9.251932720091893e-05) nohead lc rgb 'red' lw line_w
plot "spectrum.txt" index 4 u (($3 > 2.84)? ($3 - 3.23-x1) : ($3 + 3.053-x1)):(-log($2)) w p notitle pt a lc b ps c

unset arrow
unset label
set lmargin at screen 0.52
set rmargin at screen 0.6
set xlabel "$\\Delta  k_{y}$" offset kx,ky
unset ytics
set label "\\tiny $ S_{z} = 1$" at 1-x1,0.4
set label "1" at 4.0878 - d-x1,1.5 tc rgb "red"
set label "2" at 3.19 - d-x1,2.8-e tc rgb "red"
set arrow from 3.19-x1,-log(0.04207561470370) to 3.19-x1,-log(0.033058813893777) nohead lc rgb 'red' lw line_w
set label "4" at 2.293 - d-x1,3.3-e tc rgb "red"
set arrow from 2.293-x1,-log(0.0234398254893215) to 2.293-x1,-log(0.0101927995217407) nohead lc rgb 'red' lw line_w
set label "8" at 1.396 - d-x1,4.4-e tc rgb 'red'
set arrow from 1.396-x1,-log(0.008411066855) to 1.396-x1,-log(0.0023518969601) nohead lc rgb 'red' lw line_w
set label "14" at 0.498 - 2*d-x1,5.5-e tc rgb 'red'
set arrow from 0.498-x1,-log(0.002944654050485286) to 0.498-x1,-log(9.251932720091893e-05) nohead lc rgb 'red' lw line_w
plot "spectrum.txt" index 5 u (($3 > 2.84)? ($3 - 3.23-x1) : ($3 + 3.053-x1)):(-log($2)) w p notitle pt a lc b ps c

unset arrow
unset label
set lmargin at screen 0.6
set rmargin at screen 0.68
set xlabel "$\\Delta  k_{y}$" offset kx,ky
unset ytics
set label "\\tiny $ S_{z} = 2$" at 1-x1,0.4
set label "1" at 3.19 - d-x1,3 tc rgb "red"
set label "2" at 2.293 - d-x1,3.7-e tc rgb "red"
set arrow from 2.293-x1,-log(0.0169561123386) to 2.293-x1,-log(0.01019284373447120) nohead lc rgb 'red' lw line_w
set label "4" at 1.396 - d-x1,4.4-e tc rgb 'red'
set arrow from 1.396-x1,-log(0.008411066855) to 1.396-x1,-log(0.0023518969601) nohead lc rgb 'red' lw line_w
set label "8" at 0.498 - d-x1,5.5-e tc rgb 'red'
set arrow from 0.498-x1,-log(0.00294471015986) to 0.498-x1,-log(0.000364332075696) nohead lc rgb 'red' lw line_w
plot "spectrum.txt" index 6 u (($3 > 2.84)? ($3 - 3.23-x1) : ($3 + 3.053-x1)):(-log($2)) w p notitle pt a lc b ps c


#3
x2 = 4.132
reset
set ylabel "$ -ln(\\lambda _{i})$" offset 1.7
set yr [0:10.4]
set xr [-0.5-x2:6.5-x2]
set xlabel "$\\Delta  k_{y}$" offset kx,ky
set label "(c) Fermion sector" at -0.5-x2,10.65
set lmargin at screen 0.75
set rmargin at screen 0.83
set tmargin at screen 0.95
set bmargin at screen 0.1
set border lw border_w
set xtics -4,2,0 offset 0,0.5
set ytics offset 0.8
set label "\\tiny $ S_{z} = -2$" at 0.5-x2,0.4
set label "1" at 4.3465 - d-x2,1.7 tc rgb "red"
set label "1" at 3.46 - d-x2,3-e tc rgb "red"
set label "3" at 2.6746 - d-x2,3.6-e tc rgb "red"
set arrow from 2.6746-x2,-log(0.02006956246) to 2.6746-x2,-log(0.01035276545558) nohead lc rgb 'red' lw line_w
set label "5" at 1.8896 - d-x2,4.1-e tc rgb 'red'
set arrow from 1.8896-x2,-log(0.011268879408971) to 1.8896-x2,-log(0.002993454058443) nohead lc rgb 'red' lw line_w
set label "10" at 1.1245 - 2*d-x2,5.3-e tc rgb 'red'
set arrow from 1.1245-x2,-log(0.0035461341139187) to 1.1245-x2,-log(0.0005200228524734) nohead lc rgb 'red' lw line_w
set label "16" at 0.41 - 2*d-x2,5.8-e tc rgb 'red'
set arrow from 0.32-x2,-log(0.0021663664612694) to 0.32-x2,-log(4.2620256633385436e-05) nohead lc rgb 'red' lw line_w
plot "spectrum.txt" index 7 u (($3 > -2.52)? ($3 + 2.1-x2) : ($3 + 8.383-x2)):(-log($2)) w p notitle pt a lc b ps c

unset arrow
unset label
unset ylabel
set lmargin at screen 0.83
set rmargin at screen 0.91
set xlabel "$\\Delta  k_{y}$" offset kx,ky
unset ytics
set label "\\tiny $ S_{z} = -1$" at 0.45-x2,0.4
set xtics -4,2,0 offset 0,0.5
set label "1" at 4.2319 - d-x2,1.7 tc rgb "red"
set label "2" at 3.3465 - d-x2,3-e tc rgb "red"
set arrow from 3.3465-x2,-log(0.035466214672595) to 3.3465-x2,-log(0.03407282558676824) nohead lc rgb 'red' lw line_w
set label "4" at 2.56 - d-x2,3.6-e tc rgb "red"
set arrow from 2.56-x2,-log(0.02006956246) to 2.56-x2,-log(0.01035276545558) nohead lc rgb 'red'
set label "7" at 1.775 - d-x2,4.1-e tc rgb 'red'
set arrow from 1.775-x2,-log(0.011268879408971) to 1.775-x2,-log(0.002993454058443) nohead lc rgb 'red' lw line_w
set label "13" at 1.009 - 2*d-x2,5.3-e tc rgb 'red'
set arrow from 1.009-x2,-log(0.0035461341139187) to 1.009-x2,-log(0.0005200228524734) nohead lc rgb 'red' lw line_w
set label "21" at 0.206 - 2*d-x2,5.8-e tc rgb 'red'
set arrow from 0.206-x2,-log(0.0021663664612694) to 0.206-x2,-log(4.2633877514559656e-05) nohead lc rgb 'red' lw line_w
plot "spectrum.txt" index 8 u (($3 > -1)? ($3 + 1.2-x2) : ($3 + 7.483-x2)):(-log($2)) w p notitle pt a lc b ps c

unset arrow
unset label
set lmargin at screen 0.91
set rmargin at screen 0.99
set xlabel "$\\Delta  k_{y}$" offset kx,ky
unset ytics
set label "\\tiny $ S_{z} = 0$" at 1-x2,0.4
set label "1" at 4.3173 - d-x2,1.7 tc rgb "red"
set label "1" at 3.43 - d-x2,3-e tc rgb "red"
set label "3" at 2.6446 - d-x2,3.6-e tc rgb "red"
set arrow from 2.6446-x2,-log(0.02006956246) to 2.6446-x2,-log(0.01035276545558) nohead lc rgb 'red' lw line_w
set label "5" at 1.8596 - d-x2,4.1-e tc rgb 'red'
set arrow from 1.8596-x2,-log(0.011268879408971) to 1.8596-x2,-log(0.002993454058443) nohead lc rgb 'red' lw line_w
set label "10" at 1.0945 - 2*d-x2,5.3-e tc rgb 'red'
set arrow from 1.0945-x2,-log(0.0035461341139187) to 1.0945-x2,-log(0.0005200228524734) nohead lc rgb 'red' lw line_w
set label "16" at 0.38 - 2*d + 0.03-x2,5.8-e tc rgb 'red'
set arrow from 0.29-x2,-log(0.0021663664612694) to 0.29-x2,-log(4.2620256633385436e-05) nohead lc rgb 'red' lw line_w
plot "spectrum.txt" index 9 u (($3 > -0.8)? ($3 + 0.5-x2) : ($3 + 6.783-x2)):(-log($2)) w p notitle pt a lc b ps c



unset multiplot
