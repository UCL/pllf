/*
pllf_examples_2007paper.do
Reproduce all graphs in the 2007 paper
IW 3nov2025 based on files
	asym.do
	bmft_samplesize.do
	examples.do
*/

clear
set more off
set scheme sj
version 12

/*
	Example 1 
*/

* Figure 1
use brcancer, clear
stset rectime censrec
fracgen x1 -2 -0.5
fracgen x6 0.5
pllf, profile(x4a) gropt(name(fig1, replace)): stcox x1_1 x1_2 x4a x5e x6_1 hormon



* Section 5.2: effect of smaple size
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
	qui gen `v'=.
}
local i 0
qui forval j=5(5)100 {	// % of the observations
	local ++i
	local n=int((`j'+.001)*_N/100)
	sort order
	pllf, profile(x4a) nograph: stcox x1_1 x1_2 x4a x5e x6_1 hormon in 1/`n'
	foreach v in n_llci n_ulci l_llci l_ulci {
		replace `v'=r(`v') in `i'
	}
	replace p=`j' in `i'
	replace b=r(b) in `i'
	replace asym=r(asym) in `i'
	noi di `j', _c
}
keep order-asym
keep if !mi(p)

// Figure 2
line n_llci l_llci n_ulci l_ulci b p, sort clp(l - l - l) lcolor(gs8 = = = gs4) ///
 leg(label(1 "Normal-based") label(2 "PLL-based") label(5 "MLE of beta") ///
 order(1 2 5) ring(0) col(1) pos(1)) xtitle("Percentage of observations included", size(large)) ///
 ytitle("Confidence limits", size(large)) name(fig2, replace)

// Figure 3
line asym p, sort yline(0) yscale(r(0 25)) yla(0(5)25) ///
 xtitle("Percentage of observations included", size(large)) ///
 ytitle("Asymmetry (A)", size(large)) name(fig3, replace)



/*
	Example 2
*/
use brcancer, clear

pllf, formula(exp(-@*x5)) range(-1 1) gropt(name(fig4,replace)): stcox

pllf, formula(exp(-@*x5)) range(.05 .25) gropt(name(fig5,replace)): stcox

// Score function
dydx _pll _beta, gen(d)
line d _beta, sort xline(.1174) yline(0) ytitle("Score function") ///
 name(fig6, replace)
