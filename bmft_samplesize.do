set scheme s1mono
global stub bmft_samplesize
/*
	Vary the sample size and look at the asymmetry of PLL CI for x4a.
	Univariate Cox model.
*/
use \a\c233\brcancer, clear
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

save bmft_samplesize, replace
