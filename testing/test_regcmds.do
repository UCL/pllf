/*
test_regcmds.do
Test all the regression commands listed
by running the first example from each help file
IW 3nov2025
*/

/* clogit */

webuse lowbirth2
clogit low lwt smoke ptd ht ui i.race, group(pairid)
xi: pllf, profile(lwt): clogit low lwt smoke ptd ht ui i.race, group(pairid)

/* cnreg -> intreg */
webuse intregxmpl, clear
intreg wage1 wage2 age c.age#c.age nev_mar rural school tenure
gen agesq=age^2
pllf, profile(age): intreg wage1 wage2 age agesq nev_mar rural school tenure

/* glm, logistic, logit, probit */
webuse lbw, clear
glm low age lwt i.race smoke ptl ht ui, family(binomial)
xi: pllf, profile(age): glm low age lwt i.race smoke ptl ht ui, family(binomial)
xi: pllf, profile(age): logistic low age lwt i.race smoke ptl ht ui
xi: pllf, profile(age): probit low age lwt i.race smoke ptl ht ui

/* heckman */
webuse womenwk
heckman wage educ age, select(married children educ age)
pllf, profile(educ): heckman wage educ age, select(married children educ age)

/* mlogit */
webuse sysdsn1
mlogit insure age male nonwhite i.site
xi: pllf, profile([Prepaid]age): mlogit insure age male nonwhite i.site

/* nbreg, gnbreg */
webuse rod93, clear
generate logexp=ln(exposure)
nbreg deaths i.cohort, exposure(exp)
xi: pllf, profile(_Icohort_2): nbreg deaths i.cohort, exposure(exp)
gnbreg deaths age_mos, lnalpha(i.cohort) offset(logexp)
xi: pllf, verbose profile(age_mos): gnbreg deaths age_mos, lnalpha(i.cohort) offset(logexp)

/* ologit, oprobit */
webuse fullauto, clear
ologit rep77 foreign length mpg
pllf, profile(foreign): ologit rep77 foreign length mpg
pllf, profile(foreign): oprobit rep77 foreign length mpg

/* poisson */
webuse dollhill3, clear
poisson deaths smokes i.agecat, exposure(pyears)
xi: pllf, profile(smokes): poisson deaths smokes i.agecat, exposure(pyears)

/* regress */
sysuse auto, clear
regress mpg weight foreign
pllf, profile(weight): regress mpg weight foreign

/* reg3: NOT SUPPORTED */
/*
webuse klein, clear
reg3 (consump wagepriv wagegovt) (wagepriv consump govt capital1)
pllf, profile(wagepriv): reg3 (consump wagepriv wagegovt) (wagepriv consump govt capital1)
*/

/* stcox, streg, stpm, stpm2 */
webuse brcancer, clear
stset rectime, failure(censrec = 1)
stcox hormon, nohr
pllf, profile(hormon): stcox hormon, nohr
streg hormon, dist(weibull) nohr
pllf, profile(hormon): streg hormon, dist(weibull) nohr
stpm hormon, scale(hazard) df(4) 
pllf, profile(hormon): stpm hormon, scale(hazard) df(4) 
stpm2 hormon, scale(hazard) df(4)
pllf, profile(hormon): stpm2 hormon, scale(hazard) df(4) 

/* stpm3: not supported by pllf at present
since stpm3 doesn't support offset() */
* pllf, profile(hormon): stpm3 hormon, scale(lnhazard) df(4) 

/* mixed: FAILS */
/*
webuse nlswork, clear
mixed ln_w grade age c.age#c.age ttl_exp tenure c.tenure#c.tenure || id:
gen agesq=age^2
gen tenuresq=tenure^2
pllf, profile(grade): mixed ln_w grade age agesq ttl_exp tenure tenuresq || id:
*/
