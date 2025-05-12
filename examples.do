clear
set more off
set scheme sj
/*
	Example 1 for pllf paper.
*/
sjlog using pllf1, replace
use brcancer
stset rectime censrec
fracgen x1 -2 -0.5
fracgen x6 0.5
pllf stcox x1_1 x1_2 x4a x5e x6_1 hormon, profile(x4a) gropt(saving(fig1, replace))
sjlog close, replace

graph export pllf1.eps, replace


/* sjlog using pllf2, replace
global stub bmft_samplesize
/*
	Vary the sample size and look at the asymmetry of PLL CI for x4a.
	Univariate Cox model.
*/
use brcancer, clear
stset rectime censr
fracgen x1 -2 -.5
fracgen x6 .5

* Random permutation of the observations
set seed 10101
gen u=uniform()

sort u
gen int order=_n
foreach v in n_llci n_ulci l_llci l_ulci p b asym {
	gen `v'=.
}
local i 0
qui forval j=5(5)100 {	// % of the observations
	local ++i
	local n=int((`j'+.001)*_N/100)
	sort order
	pllf stcox x1_1 x1_2 x4a x5e x6_1 hormon in 1/`n', profile(x4a) nograph
	foreach v in n_llci n_ulci l_llci l_ulci {
		replace `v'=r(`v') in `i'
	}
	replace p=`j' in `i'
	replace b=r(b) in `i'
	replace asym=r(asym) in `i'
	noi di `j', _c
}

save bmft_samplesize

sjlog close, replace */


sjlog using pllf3, replace

global stub bmft_samplesize_plot
/*
	Plot results from bmft_samplesize
*/
use bmft_samplesize, clear

line n_llci l_llci n_ulci l_ulci b p, sort clp(l - l - l) lcolor(gs8 = = = gs4) ///
 leg(label(1 "Normal-based") label(2 "PLL-based") label(5 "MLE of beta") ///
 order(1 2 5) ring(0) col(1) pos(1)) xtitle("Percentage of observations included", size(large)) ///
 ytitle("Confidence limits", size(large)) sav(${stub}_1, replace)

line asym p, sort yline(0) yscale(r(0 25)) yla(0(5)25) ///
 xtitle("Percentage of observations included", size(large)) ///
 ytitle("Asymmetry (A)", size(large)) sav(${stub}_2, replace)

sjlog close, replace

graph export ${stub}_1.eps, replace
graph export ${stub}_2.eps, replace


sjlog using pllf4, replace

global stub figz
/*
	For example 2, using syntax 2
*/
use brcancer

pllf stcox X, formula(exp(-X*x5)) range(-1 1) gropt(sav(${stub}_1,replace))
return list

pllf stcox X, formula(exp(-X*x5)) range(.05 .25) gropt(sav(${stub}_2,replace))
return list

// Score function
dydx _pll _beta, gen(d)
line d _beta, sort xline(.1174) yline(0) ytitle("Score function") ///
 sav(${stub}_3, replace)

sjlog close, replace

graph export ${stub}_1.eps, replace
graph export ${stub}_2.eps, replace
graph export ${stub}_3.eps, replace

