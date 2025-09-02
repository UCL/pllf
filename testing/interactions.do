/*
Can we profile an interaction?
Yes, but only using xi?
*/

pda
use ../brcancer, clear
stcox x2##x4a x5e x6 hormon, nohr

// Can we profile the interaction parameter?

* not using new factor variables
* but we can using xi

* the easy way
xi i.x2*i.x4a
pllf, profile(_Ix2Xx4a_2_1) normal n_eval(20): stcox _I* x5e x6 hormon, nohr

* the one-step way
xi: pllf, profile(_Ix2Xx4a_2_1) normal n_eval(20): stcox i.x2*i.x4a x5e x6 hormon