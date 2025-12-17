/*
Test that the code handles collinearity
test_collinear.do
IW 23jul2025
*/

prog drop _all

clear
input Z Y n
0 0 38
0 1 4
1 0 36
1 1 1
end
tab Z Y [fw=n]

* correct command
pllf, profile(Z): logit Y Z [fw=n]
local ref = r(l_llci)

* duplicated covariate: fails by default
cap noi pllf, profile(Z): logit Y Z Z [fw=n]
assert _rc==498

* duplicated covariate: succeeds with dropcol option
pllf, profile(Z) dropcol: logit Y Z Z [fw=n]
assert `ref' == r(l_llci)

* identical covariate: fails by default
gen Z2=Z
cap noi pllf, profile(Z): logit Y Z2 Z [fw=n]
assert _rc==498

* identical covariate placed first: fails with dropcol option
cap noi pllf, profile(Z) dropcol: logit Y Z2 Z [fw=n]
assert _rc==198

* identical covariate placed last: succeeds with dropcol option
pllf, profile(Z) dropcol: logit Y Z Z2 [fw=n]
assert `ref' == r(l_llci)
