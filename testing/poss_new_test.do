/*
Can we profile an interaction?
Yes, but only using xi?
*/

pda
use ../brcancer, clear
stcox x2##x4a x5e x6 hormon, nohr
xi i.x2*i.x4a
pllf, profile(_Ix2Xx4a_2_1) shownormal n_eval(20): stcox _I* x5e x6 hormon, nohr
