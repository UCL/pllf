/*
Main test file for pllf
IW 17mar2025
Added commands from help file 12may2025
Revised for the colon command 28may2025
Added tests of normal() and collinearity handling 1aug2025
Update for renaming shownormal to normal, 1sep2025
Test all graph options, 10dec2025
*/

* from help file
use ../brcancer, clear

* Syntax 1

stcox x1 x4a x5e x6 hormon, nohr
local b5 = _b[x5e]
local se5 = _se[x5e]
pllf, profile(x5e) range(-3 -1) normal: stcox x1 x4a x5e x6 hormon, nohr
assert `b5' == r(b)
assert `se5' == r(se)
assert reldif(r(se),r(pse))<0.001

stpm x1 x4a x5e x6 hormon, df(2) scale(h) 
local b1 = [xb]_b[x1]
local se1 = [xb]_se[x1]
pllf, profile(x1) gen(betavar pllvar) normal: stpm x1 x4a x5e x6 hormon, df(2) scale(h) 
confirm var betavar pllvar _pllnorm
mac l _b1
di r(b)
assert reldif(`b1',r(b))<1E-7
assert reldif(`se1',r(se))<1E-7
assert reldif(r(se),r(pse))<0.001
drop betavar pllvar _pllnorm

* with/without eqname: failed 1/8/2025
* note it's faster via constraints
pllf, profile(x5e) n_eval(20): streg x1 x4a x5e x6 hormon, dist(expo)
local ulci = r(l_ulci)
pllf, profile([_t]x5e) n_eval(20): streg x1 x4a x5e x6 hormon, dist(expo)
assert reldif(`ulci', r(l_ulci))<1E-7

* collinearity
gen horm2=hormon
cap noi pllf, profile(x5e) n_eval(20): streg horm2 x1 x4a x5e x6 hormon, distribution(expo)
assert _rc==498
pllf, profile(x5e) n_eval(20) dropcol: streg horm2 x1 x4a x5e x6 hormon, distribution(expo)
assert reldif(`ulci', r(l_ulci))<1E-7

pllf, profile([ln_p]_cons) n_eval(50): streg x1 x4a x5e x6 hormon, distribution(weibull)

pllf, profile([ln_p]x4b) deviance difference n_eval(20): streg x1 x4a x5e x6 hormon, distribution(weibull) ancillary(x4b) 


* Syntax 2
pllf, formula(exp(-@*x5)) range(.05 .25): stcox x1 x4a x6 hormon
local pllf_ll = r(ll)
mac l _pllf_ll
gen expXx5 = exp(-r(b)*x5)
stcox x1 x4a expXx5 x6 hormon
di e(ll)
assert abs(e(ll)-`pllf_ll')<1E-3
pllf, formula(exp(-#*x5)) placeholder(#) range(.05 .25): stcox x1 x4a x6 hormon
di r(ll)
assert abs(r(ll)-`pllf_ll')<1E-3

* and the case where we need to differentiate
pllf, trace formula(exp(-@*x5/100)) n(5) range(-20 20): stcox x1 x5 x4a x6 hormon

* Syntax 1 again
sysuse auto, clear
gen one = 1
pllf, profile(one): reg mpg foreign rep78 one, noconstant 
local asym=r(asym)
pllf, profile(_cons): reg mpg foreign rep78
assert reldif( `asym', r(asym)) < 1E-7


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
pllf, profile(z) range(-5 0) eform: logit d z [fw=n]


* blogit and glm
clear
input z n d
0 100 38
1 100 28
end

blogit d n i.z
pllf, profile(z) debug gropt(name(blogit,replace)) normal(lcol(red)): blogit d n z
* understands blogit has 2 "yvars"

glm d z, family(binomial n) 
pllf, profile(z) debug: glm d z, family(binomial n)

* various graph options
pllf, profile(z) mleline: glm d z, family(binomial n)
pllf, profile(z) cilines(lcol(red)) levline(lcol(blue)) mleline: glm d z, family(binomial n)
pllf, profile(z) cilines(off) levline(off): glm d z, family(binomial n)
pllf, profile(z) nograph nodots: glm d z, family(binomial n)

* poisson with exposure()
clear
input group pyears events
0 200 38
1 100 19
end
poisson events group, exposure(pyears)
pllf, profile(group): poisson events group, exposure(pyears)
pllf, profile(group) range(-1 1): poisson events group, exposure(pyears)
local asym=r(asym)
di r(asym)
pllf, profile([events]group) range(-1 1): poisson events group, exposure(pyears)
di r(asym)
assert reldif( `asym', r(asym)) < 1E-4 // fails at 1E-7

* profile the constant: two ways
pllf, trace n_eval(10) debug profile(_cons): poisson events group, exposure(pyears)
local asym=r(asym)
gen one=1
pllf, trace n_eval(10) debug profile(one): poisson events group one, exposure(pyears) nocons
assert reldif( `asym', r(asym)) < 1E-7


* TRISST data, MRI vs CT comparison
* NB it's a NI trial so don't get excited about the profile CI not crossing the null 
* 	NI margin is +0.057
use ../TRISST, clear
pllf, norm verbose profile(modality): binreg outcome modality [fw=n], rd


* TEST HANDLING INTERACTIONS
use ../brcancer, clear
stcox x2##x4a x5e x6 hormon, nohr

* we can't profile the interaction parameter using new factor variables,
* but we can using xi

* the easy way
xi i.x2*i.x4a
pllf, profile(_Ix2Xx4a_2_1) normal n_eval(20): stcox _I* x5e x6 hormon, nohr

* the one-step way
xi: pllf, profile(_Ix2Xx4a_2_1) normal n_eval(20): stcox i.x2*i.x4a x5e x6 hormon
