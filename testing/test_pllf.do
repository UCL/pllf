/*
rudimentary test file for pllf
IW 17mar2025
Added commands from help file 12may2025
REvised for the colon command 28may2025
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

stcox x1 x4a x5e x6 hormon, nohr
local b5 = _b[x5e]
local se5 = _se[x5e]
pllf, profile(x5e) range(-3 -1): stcox x1 x4a x5e x6 hormon, nohr
assert `b5' == r(b)
assert `se5' == r(se)
assert reldif(r(se),r(pse))<0.001

stpm x1 x4a x5e x6 hormon, df(2) scale(h) 
local b1 = [xb]_b[x1]
local se1 = [xb]_se[x1]
pllf, profile(x1) gen(X Y) normcoll: stpm x1 x4a x5e x6 hormon, df(2) scale(h) 
confirm var X Y
mac l _b1
di r(b)
assert reldif(`b1',r(b))<1E-7
assert reldif(`se1',r(se))<1E-7
assert reldif(r(se),r(pse))<0.001

pllf, profile([ln_p]_cons) n(50): streg x1 x4a x5e x6 hormon, distribution(weibull)

pllf, profile([ln_p]x4b) deviance difference n(20): streg x1 x4a x5e x6 hormon, distribution(weibull) ancillary(x4b) 

* pllf, profile([d]group) range(0.27 1.55): poisson d group, exposure(y) 

* Syntax 2
*pllf logit y x1 X, formula(exp(-X*x2)) range(.05 .25)
pllf, formula(exp(-X*x5)) range(.05 .25): stcox x1 x4a X x6 hormon
local pllf_ll = r(ll)
mac l _pllf_ll
gen expXx5 = exp(-r(b)*x5)
stcox x1 x4a expXx5 x6 hormon
di e(ll)
assert abs(e(ll)-`pllf_ll')<1E-3

* Syntax 1 again
sysuse auto, clear
gen one = 1
pllf, profile(one) normcoll: logit foreign mpg one, noconstant 
* NB normcoll is needed here


* perfect prediction
clear
input z d n
0 0 38
0 1 4
1 0 36
1 1 0
end
tab d z [fw=n], chi2
pllf, profile(z) range(-5 0) gropt(name(perfpred,replace)):  logit d z [fw=n]



* blogit and glm
clear
input z n d
0 100 38
1 100 28
end

blogit d n i.z
pllf, profile(z) debug normcoll gropt(name(blogit,replace)):  blogit d n z
* understands blogit has 2 "yvars"

glm d z, family(binomial n) 
pllf, profile(z) debug normcoll: glm d z, family(binomial n)

di as result "*** PLLF HAS PASSED ALL THE TESTS IN `filename'.do ***"

log close
