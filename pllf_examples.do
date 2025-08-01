clear
input group pyears events
      0 200 38
      1 100 19
      end

* Fit Poisson model

poisson events group, exposure(pyears)

* Explore profile likelihood for the coefficient of group

pllf, profile(group): poisson events group, exposure(pyears)

* Load TRISST trial data (Joffe et al, 2022). The data 
* are from a non-inferiority trial of MRI vs CT for surveillance after testicular 
* cancer. The PLL CI does not cross zero while the Normal CI does. However this 
* is a non-inferiority trial with margin a risk difference of +0.057, so both CIs 
* clearly establish non-inferiority.

use TRISST, clear
pllf, shownorm verbose profile(modality): binreg outcome modality [fw=n], rd

* Load breast cancer data

webuse brcancer, clear
stset rectime, failure(censrec) scale(365.24)

* Explore profile likelihood for coefficient of x5e

pllf, profile(x5e) range(-3 -1): stcox x1 x4a x5e x6 hormon

* Explore profile likelihood for coefficient of x1

pllf, profile(x1) gen(X Y): stpm2 x1 x4a x5e x6 hormon, df(2) scale(hazard) 

* Explore profile likelihood for the constant

pllf, profile(_cons) n(50): streg x1 x4a x5e x6 hormon, distribution(weibull)

* Explore profile likelihood for Weibull shape parameter

pllf, profile([ln_p]_cons) n(50): streg x1 x4a x5e x6 hormon, distribution(weibull)

* Explore profile likelihood for predictor of Weibull shape parameter

pllf, profile([ln_p]x4b) deviance difference n(20): streg x1 x4a x5e x6 hormon, distribution(weibull) ancillary(x4b) 

* Syntax 2

* The following two commands are equivalent.

pllf, formula(exp(-X*x5)) range(.05 .25): stcox x1 x4a X x6 hormon
pllf, placeholder(@) formula(exp(-@*x5)) range(.05 .25): stcox x1 x4a @ x6 hormon
