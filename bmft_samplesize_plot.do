set scheme s1mono
global stub bmft_samplesize_plot
/*
	Plot results from bmft_samplesize
*/
use bmft_samplesize, replace

line n_llci l_llci n_ulci l_ulci b p, sort clp(l - l - l) lcolor(gs8 = = = gs4) ///
 leg(label(1 "Normal-based") label(2 "PLL-based") label(5 "MLE of beta") ///
 order(1 2 5) ring(0) col(1) pos(1)) xtitle("Percentage of observations included", size(large)) ///
 ytitle("Confidence limits", size(large)) sav(${stub}_1, replace)

line asym p, sort yline(0) yscale(r(0 25)) yla(0(5)25) ///
 xtitle("Percentage of observations included", size(large)) ///
 ytitle("Asymmetry (A)", size(large)) sav(${stub}_2, replace)

gph2emf ${stub}_1
gph2emf ${stub}_2
