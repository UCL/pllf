
myadopath pllf
pda

webuse nlswork, clear
*mixed ln_w grade age c.age#c.age ttl_exp tenure c.tenure#c.tenure || id:
gen agesq=age^2
gen tenuresq=tenure^2

// this doesn't work, because mixed doesn't allow constraints...
pllf, verbose trace profile([ln_wage]grade) normcoll: mixed ln_w grade age agesq ttl_exp tenure tenuresq || id:

// ... help file says this shoudl work, but it doesn't
constraint 3 [ln_wage]grade=.0665654765887296
mixed ln_w grade age agesq ttl_exp tenure tenuresq, constraints(3) || id:

// this works!
pllf, verbose n(11) debug profile([ln_wage]grade) normcoll: meglm ln_w grade age agesq ttl_exp tenure tenuresq || id:

