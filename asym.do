global stub asym
/*
	Asymmetry of each predictor.
*/
use brcancer, replace
stset rectime censrec
fracgen x1 -2 -0.5
fracgen x6 0.5
local vars x1_1 x1_2 x4a x5e x6_1 hormon
foreach v in `vars' {
	qui pllf stcox `vars', profile(`v') nograph
	di "`v'", r(asym)
}
