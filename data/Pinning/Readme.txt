spin-1 J1-J2-Jchi square model main paper Fig. 1 (right figure)

The file "Pining_diff_phase_hz01.txt" contains the onsite <S^z_i> when applying a edge pinning field hz=0.1 for four different parameter choices in four phases labeled in the file. 
The first column stands for the site i, and the second column the <S^z_i>

The files "KspaceSpinStru_Ny_X_Ny_Nx24Ny8j2{J2}jz1jk{Jchi}omega0SU2.txt" contains the spin structure factor for Jchi=Jchi, J2=J2 summing over the middle Ly*Ly sites. The results are obtained with SU(2) codes on 24*8 cylinders. The definition of the spin structure factor can be found in the caption of Fig. 1
The first column stands for kx, the second column ky, and the third column S(kx,ky)


To compile, use
$gnuplot Gnuplot_Pining_four_phase.sh
$latex Pining_four_phase.tex
$dvips Pining_four_phase.dvi
$ps2pdf Pining_four_phase.ps
