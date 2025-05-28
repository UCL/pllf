/*
rudimentary test file for pllf
IW 17mar2025
Added commands from help file 12may2025
REvised for the colon command 28may2025
*/

local filename test_pllf_colon

prog drop _all
cd c:\ian\git\pllf\testing
cap log close
set linesize 100
clear all // avoids the "too many sersets" error

log using `filename', replace text
version
which pllf

// START TESTING

* from help file
use ../brcancer, clear
* Syntax 1
pllf, profile(x5e) range(-3 -1): stcox x1 x4a x5e x6 hormon
pllf, profile(x1): stpm x1 x4a x5e x6 hormon, df(2) scale(h) gen(X Y) 
pllf, profile([ln_p]_cons) n(50): streg x1 x4a x5e x6 hormon, distribution(weibull)
pllf, profile([ln_p]x4b) deviance difference n(20): streg x1 x4a x5e x6 hormon, distribution(weibull) ancillary(x4b) 
* pllf, profile([d]group) range(0.27 1.55): poisson d group, exposure(y) 

* Syntax 2
*pllf logit y x1 X, formula(exp(-X*x2)) range(.05 .25)
pllf, formula(exp(-X*x5)) range(.05 .25): stcox x1 x4a X x6 hormon

* Syntax 1 again
sysuse auto, clear
gen one = 1
pllf, profile(one): logit foreign mpg one, noconstant 



* perfect prediction
clear
input z d n
0 0 38
0 1 4
1 0 36
1 1 0
end
tab d z [fw=n]
pllf, profile(z) range(-5 0):  logit d z [fw=n]



* blogit and glm
clear
input z n d
0 100 38
1 100 28
end

blogit d n i.z
pllf, profile(z) debug normcoll:  blogit d n z
* doesn't understand blogit has 2 "yvars"

glm d z, family(binomial n) 
pllf, profile(z) debug normcoll: glm d z, family(binomial n)

di as result "*** PLLF HAS PASSED ALL THE TESTS IN `filename'.do ***"

log close
