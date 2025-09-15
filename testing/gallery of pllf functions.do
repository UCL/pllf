/* 
Gallery of PLLF functions
IW 12sep2025
*/
set seed 4561964
clear

// t-test
set obs 10
gen z = _n>_N/2
gen y = 0.5*z + rnormal()
reg y z
pllf, profile(z) norm range(-5 5) gropt(name(reg,replace)): reg y z

// logistic
gen yb = inrange(_n,5,8)
pllf, profile(z) norm range(-3 7) gropt(name(logit,replace)): logit yb z

// Cox
gen t = -log(runiform())
stset t
pllf, profile(z) norm range(-5 8) gropt(name(cox,replace)): stcox z

// Complex
clear
set obs 100
drawnorm x1-x10
gen t = -5*log(runiform())/exp(0.3*(x1+x2-x3-x4))
gen d = t>=5
replace t=1 if !d
tab d
stset t, fail(d)
pllf, profile(x1) norm gropt(name(complex,replace)): stcox x1-x10
