set scheme s1mono
global stub figz
/*
	For example 2, using syntax 2
*/
set more off
logopen
use brcancer, clear

pllf stcox X, formula(exp(-X*x5)) range(-1 1) gropt(sav(${stub}_1,replace))
return list

pllf stcox X, formula(exp(-X*x5)) range(.05 .25) gropt(sav(${stub}_2,replace))
return list

// Score function
dydx _pll _beta, gen(d)
line d _beta, sort xline(.1174) yline(0) ytitle("Score function") ///
 sav(${stub}_3, replace)

gph2emf ${stub}_1
gph2emf ${stub}_2
gph2emf ${stub}_3

log close
