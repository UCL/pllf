/*
Explore PLLFs
explore_pllf.do
3 data sets showing RR=0.5 and same Wald CI
IW 8nov2024

I could use this with a 2nd study with symmetrical LL, e.g. 3/10 vs 3/10.
I could show these 3 data sets give the same result when 2-stage analysed using normal approximation, but not when correctly analysed
*/
pda

clear
input z d1 py1 d2 py2 d3 py3
0 3 10 2 10 6 30
1 3 20 6 60 2 20
end

forvalues i=1/3 {
	di as input "Data data`i'"
	gen lpy`i' = log(py`i')
	poisson d`i' z, offset(lpy`i') irr
	local b`i'=_b[z]
	local se`i'=_se[z]
	est store d`i'
	pllf, profile(z) range(-3 1.5) gen(beta`i' pll`i') nograph difference: poisson d`i' z, offset(lpy`i') irr 
	local d0 = d`i'[1]
	local d1 = d`i'[2]
	local py0 = py`i'[1]
	local py1 = py`i'[2]
	label var pll`i' "Data `i': `d1'/`py1' vs `d0'/`py0'"
}

* Verify Wald results are the same
est table _all, eq(1) keep(z) se

* Create Normal approximation to LL
* same for all data sets
gen normapprox_ll = -0.5*(beta1-`b1')^2/(`se1')^2
label var norma "Normal approx"

* Graph
line pll* norma beta1, ytitle(Shifted PLLF) xtitle(log RR) yli(-1.92) saving(explore_pllf, replace) lcol(red blue green black)
