/*
rudimentary test file for pllf
IW 17mar2025
Added commands from help file 12may2025
*/

local filename test_pllf

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
pllf stcox x1 x4a x5e x6 hormon, profile(x5e) range(-3 -1)
pllf stpm x1 x4a x5e x6 hormon, df(2) scale(h) gen(X Y) profile(x1)
pllf streg x1 x4a x5e x6 hormon, distribution(weibull) profile([ln_p]_cons) n(50)
pllf streg x1 x4a x5e x6 hormon, distribution(weibull) ancillary(x4b) profile([ln_p]x4b) deviance difference n(20)
* pllf poisson d group, exposure(y) profile([d]group) range(0.27 1.55)

* Syntax 2
*pllf logit y x1 X, formula(exp(-X*x2)) range(.05 .25)
pllf stcox x1 x4a X x6 hormon, formula(exp(-X*x5)) range(.05 .25)

* Syntax 1 again
sysuse auto, clear
gen one = 1
pllf logit foreign mpg one, noconstant profile(one)



* perfect prediction
clear
input z d n
0 0 38
0 1 4
1 0 36
1 1 0
end
tab d z [fw=n]
pllf logit d z [fw=n], profile(z) range(-5 0)

* blogit and glm
clear
input z n d
0 100 38
1 100 28
end

blogit d n i.z
pllf blogit d n z, profile(z) debug normcoll
* doesn't understand blogit has 2 "yvars"

glm d z, family(binomial n) 
pllf glm d z, family(binomial n) profile(z) debug normcoll

di as result "*** PLLF HAS PASSED ALL THE TESTS IN `filename'.do ***"

log close
