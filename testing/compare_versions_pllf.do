/* 
CHECK OLD AND NEW VERSIONS AGREE
compare_versions_pllf
IW 23jul2025
*/

// OLD VERSION

cap adopath - $path
which pllf
pda

webuse brcancer, clear
stset rectime, fail(censrec)
stcox x1 x4a x5e x6 hormon, nohr
pllf stcox x1 x4a x5e x6 hormon, nohr profile(x5e) range(-3 -1) gen(beta1 pll1) gropt(name(pll1, replace))


// NEW VERSION

adopath ++ $path
which pllf
pda

pllf, profile(x5e) range(-3 -1) mleline gen(beta2 pll2) gropt(name(pll2, replace)): stcox x1 x4a x5e x6 hormon, nohr

assert beta1==beta2
assert abs(pll1-pll2)<1E-4 if !mi(beta1)

