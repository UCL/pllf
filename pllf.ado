*! version 1.0.1 PR 17jul2007
program define pllf, rclass sortpreserve
version 9.0
/*
	Currently supported commands include at least the following:

	clogit cnreg cox ereg fit glm gnbreg heckman logistic logit	///
		   mlogit nbreg ologit oprobit poisson probit regress reg3	///
		   streg stcox stpm weibull
*/
gettoken cmd 0 : 0
if "`cmd'"=="stpm" local eqxb [xb]
else if ("`cmd'"=="fit") | ("`cmd'"=="reg") | (substr("`cmd'",1,4)=="regr") local cmd regress

syntax anything [if] [in] [aweight fweight pweight iweight], ///
 [DEViance FORMula(string) gen(string) DIFFerence LEVel(cilevel) ///
 PLaceholder(string) PROfile(string) range(string) ///
 MAXCost(int -1) n(integer 100) noci noDOTs nograph gropt(string asis)		///
 LEVLINe(string asis) CILINes(string asis) *]

if `maxcost'<0 local maxcost = int(`n'/2)

if "`placeholder'"!="" {
	if wordcount("`placeholder'")!=1 {
		di as err "invalid placeholder()"
		exit 198
	}
}
else local placeholder X

if "`formula'"!="" {
	if "`profile'"!="" {
		di as txt "[profile() ignored]"
		local profile
	}
	if "`range'"=="" {
		di as err "range() required"
		exit 198
	}
	// Check for `placeholder' in `anything'
	local result: subinstr local anything "`placeholder'" "", count(local nat)
	if `nat'==0 {
		di as err `"`placeholder' not found in regression_cmd_stuff ( `anything' )"'
		exit 198
	}
}
else if "`profile'"!="" {
	local varlist "`anything'"
	// Disentangle profile; extract eq from it, if present
	tokenize `profile', parse("[]")
	if "`5'"!="" {
		di as err "syntax error in profile(`profile'), invalid equation name"
		exit 198
	}
	if "`2'"!="" {
		// Rudimentary check that user has entered the eq correctly
		if "`1'"!="[" | "`3'"!="]" {
			di as err "syntax error in profile(`profile'), invalid equation name"
			exit 198
		}
		local profile `4'
		local eq [`2']
		constraint free
		local cuse `r(free)'
	}
	if "`profile'"!="_cons"{
		unab profile: `profile'
	}
	else {
		if "`eq'"=="" {
			di as err "profile log likelihood for _cons not supported"
			exit 198
		}
	}
}
else {
	di as err "must supply either formula() or profile()"
	exit 198
}

if "`gen'"=="" {
	local gen1 _beta
	local gen2 _pll
}
else gettoken gen1 gen2: gen
if "`gen2'"=="" local gen2 _pll

if "`range'"!="" {
	gettoken from to: range
	confirm num `from'
	confirm num `to'
	if `from' > `to' {
		local temp `from'
		local from `to'
		local to `temp'
		local temp
	}
}

if "`profile'"!="" { // ------------ begin linear profiling --------
	marksample touse
	local wt [`weight'`exp']

	// Fit model and get level% ci. Program terminates if invalid cmd attempted.
	if "`eq'"=="" {
		qui _rmcoll `varlist' `profile'	// strips `profile' if already mentioned in `varlist'
		local tempvl `r(varlist)'
	}
	else local tempvl `varlist'
	capture noisily `cmd' `tempvl' if `touse' `wt', `options'
	local ytitle `e(depvar)'
	quietly {
		// Check that alleged parameter exists. One or both of `eqxb' or `eq' will always be null.
		capture local b0 = `eq'`eqxb'_b[`profile']
		if "`b0'"=="" local b0 .
		if "`cmd'"=="regress" local z = invttail(e(df_r), (100-`level')/200)
		else local z = -invnorm((100-`level')/200)
	
		local nobs = e(N)
		local ll0  = e(ll)
		capture local se = `eq'`eqxb'_se[`profile']
		if _rc==0 {
			local llci = `b0'-`z'*`se'
			local ulci = `b0'+`z'*`se'
		}
		else {
			local se .
			if "`range'"=="" {
				noi di as err "could not select range for `profile' (could not estimate MLE). try supplying range()"
				exit 198
			}
		}
		if "`range'"=="" {
			// previously used default range as Wald-based `level' CI; now using +/-(z*1.2)SE for this.
			local from = `b0'-`z'*1.2*`se'
			local to =   `b0'+`z'*1.2*`se'
		}

		// If no equation specified, unabbreviate varlist and remove `profile' from it
		if "`eq'"=="" {
			unab varlist: `varlist'
			local varlist: list varlist - profile
		}

		if "`cmd'"=="regress" {
			constraint free
			local cuse `r(free)'

			// For regress, need to identify yvar; otherwise, not required
			gettoken yvar varlist: varlist
		}
		if `n'>_N {
			di as txt "[Increasing the dataset size to `n'] " _cont
			set obs `n'
		}
		tempvar X Y order // values of X = regression coefficient, Y = profile likelihood for `profile'
		gen `X' = .
		gen `Y' = .
		gen long `order' = _n
		local stepsize = (`to'-`from')/(`n'-1)
		forvalues i=1/`n' {
			local b = `from'+(`i'-1)*`stepsize'
			if "`eq'"!="" {
				// Use constrained regression
				constraint define `cuse' `eq'`profile'=`b'
				`cmd' `varlist' if `touse' `wt', `options' constraint(`cuse')
			}
			else {
				if `i'==1 {
					tempvar offset
					gen `offset' = .
				}
				if "`cmd'"=="regress" {
					replace `offset' = `yvar'-`b'*`profile'
					regress `offset' `varlist' if `touse' `wt', `options'
				}
				else {
					replace `offset' = `b'*`profile'
					`cmd' `varlist' if `touse' `wt', `options' offset(`offset')
				}
			}
			sort `order'
			replace `X' = `b' in `i'
			replace `Y' = e(ll) in `i'
			if "`dots'"!="nodots" noi di "." _c
		}
		local cost 0	// "cost": number of extra evaluations of pll needed to find likelihood based CI
		local left_limit .
		local right_limit .
		if "`ci'"!="noci" {
/*
			Search for likelihood based CI bounds: VERY crude!
			Left limit first.
*/
			local target = `ll0'-`z'^2/2	// pll value for computing likelihood based CI
			if `Y'[1]<`target' {
/*
			First evaluated pll is below target pll on left of mle.
			Bracket target from already known values of pll and interpolate.
*/
				forvalues i=2/`n' {
					if `Y'[`i']>=`target' {
						local left_limit = `X'[`i'-1]+`stepsize'*(`target'-`Y'[`i'-1])/(`Y'[`i']-`Y'[`i'-1])
						continue, break	// exit forvalues loop
					}
				}
			}
			else {
/*
			Search for left ll-based confidence limit to the left of first value of b.
			Requires new evaluations of ll - no more than `maxcost' allowed.
*/
				local Yold = `Y'[1]
				local bold `from'
				local i 1
				while missing(`left_limit') & `i'<=`maxcost' {
					local b = `from'-`i'*`stepsize'
					// evaluate pll
					if "`eq'"!="" {
						// Use constrained regression
						constraint define `cuse' `eq'`profile'=`b'
						`cmd' `varlist' if `touse' `wt', `options' constraint(`cuse')
					}
					else {
						if "`cmd'"=="regress" {
							replace `offset' = `yvar'-`b'*`profile'
							regress `offset' `varlist' if `touse' `wt', `options'
						}
						else {
							replace `offset' = `b'*`profile'
							`cmd' `varlist' if `touse' `wt', `options' offset(`offset')
						}
					}
					local Ynew = e(ll)
					if `Ynew'<`target' {	// now bracketing target
						local left_limit = `bold'-`stepsize'*(`target'-`Yold')/(`Ynew'-`Yold')
					}
					else {
						local Yold `Ynew'
						local bold `b'
						local ++i
					}
					if "`dots'"!="nodots" noi di "." _c
				}
				local cost `i'
				if missing(`left_limit') noi di as txt _n "[note: failed to find left-hand confidence limit]"
			}
/*
		Now right limit
*/
			if `Y'[`n']>`target' {	// search for right_limit to the right of last value of b
				local Yold = `Y'[`n']
				local bold `to'
				local i 1
				while missing(`right_limit') & `i'<=`maxcost' {
					local b = `to'+`i'*`stepsize'
					// evaluate pll
					if "`eq'"!="" {
						// Use constrained regression
						constraint define `cuse' `eq'`profile'=`b'
						`cmd' `varlist' if `touse' `wt', `options' constraint(`cuse')
					}
					else {
						if "`cmd'"=="regress" {
							replace `offset' = `yvar'-`b'*`profile'
							regress `offset' `varlist' if `touse' `wt', `options'
						}
						else {
							replace `offset' = `b'*`profile'
							`cmd' `varlist' if `touse' `wt', `options' offset(`offset')
						}
					}
					local Ynew = e(ll)
					if `Ynew'<`target' {	// now bracketing target
						local right_limit = `bold'+`stepsize'*(`target'-`Yold')/(`Ynew'-`Yold')
					}
					else {
						local Yold `Ynew'
						local bold `b'
						local ++i
					}
					if "`dots'"!="nodots" noi di "." _c
				}
				local cost = `cost'+`i'
			}
			else {	// pll is below target on right of mle, bracket target and interpolate
				local n1 = `n'-1
				forvalues i=`n1'(-1)1 {
					if `Y'[`i']>=`target' {
						local right_limit = `X'[`i']+`stepsize'*(`target'-`Y'[`i'])/(`Y'[`i'+1]-`Y'[`i'])
						continue, break	// exit forvalues loop
					}
				}
			}
			if missing(`right_limit') noi di as txt _n "[note: failed to find right-hand confidence limit]"
		}
		lab var `X' "`eq'_b[`profile']"
	}
}
else {	// --------------- begin non-linear profiling ---------------
	marksample touse, novarlist
	local wt [`weight'`exp']
	tempvar xx
	qui gen `xx' = .
	// Trial fit of model with central value of param. Program terminates if invalid cmd attempted.
	local A = (`to'-`from')/2
	parsat `"`anything'"' if `touse', formula(`formula') var(`xx') value(`A') placeholder(`placeholder')
	cap `cmd' `r(result)' if `touse' `wt', `options'
	local rc = _rc
	if `rc' error `rc'

	local ytitle `e(depvar)'
	cap local ll = e(ll)
	if _rc | ("`ll'"==".") {
		di as err "valid log likelihood not returned in e(ll)"
		exit 198
	}

	quietly {
		if "`cmd'"=="regress" local z = invttail(e(df_r), (100-`level')/200)
		else local z = -invnorm((100-`level')/200)

		local nobs = e(N)
		replace `touse' = e(sample)

		if `n'>_N {
			di as txt "[Increasing the dataset size to `n'] " _cont
			set obs `n'
		}
		tempvar X Y order // values of X = regression coefficient, Y = profile likelihood for `profile'
		gen `X' = .
		gen `Y' = .
		gen long `order' = _n
		local stepsize = (`to'-`from')/(`n'-1)
		local y1 .
		local y2 .
		local y3 .
		local b0 .
		local ll0 .
		local done 0
		forvalues i=1/`n' {
			local A = `from'+(`i'-1)*`stepsize'
			// Substitute A in `formula' and fit model
			parsat `"`anything'"' if `touse', formula(`formula') var(`xx') value(`A') placeholder(`placeholder')
			cap `cmd' `r(result)' if `touse' `wt', `options'
			if _rc {
				noi di as err "model fit failed at parameter = " `A'
				exit 198
			}
			local y3 = e(ll)
			if !missing(`y1') & !missing(`y2') {
				if `y1'<`y2' & `y2'>`y3' {
					// Solve quadratic y = a + b*x + c*x^2 through 3 points to get MLE b0 and ll at MLE ll0
					local c = (2*(`y3'-`y2')-(`y3'-`y1'))/(2*`stepsize'^2)
					local b = (`y3'-`y2')/`stepsize'-`c'*(2*`A'-`stepsize')
					local a = `y3'-`A'*(`b'+`c'*`A')
					local b0 = -`b'/(2*`c')
					local ll0 = `a'+`b0'*(`b'+`c'*`b0')
					local done 1
				}
			}
			local y1 `y2'
			local y2 `y3'
			sort `order'
			replace `X' = `A' in `i'
			replace `Y' = e(ll) in `i'
			if "`dots'"!="nodots" noi di "." _c
		}
		if !`done' {
			// MLE not straddled. Attempt quadratic solution using terminals of range and a midpoint. Maintain equal spacing.
			if mod(`n',2)==0 {	// even number of evaluations
				local mid = `n'/2
				local last `n'-1
			}
			else {	// odd number of evaluations
				local mid = (`n'+1)/2
				local last `n'
			}
			local A = `X'[`last']
			local s = `X'[`mid']-`X'[1]	// stepsize for this exercise
			local y1 = `Y'[1]
			local y2 = `Y'[`mid']
			local y3 = `Y'[`last']
			local c = (2*(`y3'-`y2')-(`y3'-`y1'))/(2*`s'^2)
			local b = (`y3'-`y2')/`s'-`c'*(2*`A'-`s')
			local a = `y3'-`A'*(`b'+`c'*`A')
			local b0 = -`b'/(2*`c')
			local ll0 = `a'+`b0'*(`b'+`c'*`b0')
			noi di as txt _n "Note: range does not include MLE of parameter - estimate may be inaccurate."
		}
		local cost 0	// "cost": number of extra evaluations of pll needed to find likelihood based CI
		local left_limit .
		local right_limit .
		if "`ci'"!="noci" {
/*
		Search for likelihood based CI bounds: VERY crude!
		Left limit first.
*/
			local target = `ll0'-`z'^2/2	// pll value for computing likelihood based CI
			if `Y'[1]<`target' {
/*
			First evaluated pll is below target pll on left of mle.
			Bracket target from already known values of pll and interpolate.
*/
				forvalues i=2/`n' {
					if `Y'[`i']>=`target' {
						local left_limit = `X'[`i'-1]+`stepsize'*(`target'-`Y'[`i'-1])/(`Y'[`i']-`Y'[`i'-1])
						continue, break	// exit forvalues loop
					}
				}
			}
			else {
/*
			Search for left ll-based confidence limit to the left of first value of b.
			Requires new evaluations of ll - no more than `maxcost' allowed.
*/
				local Yold = `Y'[1]
				local bold `from'
				local i 1
				while missing(`left_limit') & `i'<=`maxcost' {
					local A = `from'-`i'*`stepsize'
					// evaluate pll
					parsat `"`anything'"' if `touse', formula(`formula') var(`xx') value(`A') placeholder(`placeholder')
					cap `cmd' `r(result)' if `touse' `wt', `options'
					if _rc {
						noi di as err "model fit failed at parameter = " `A'
						exit 198
					}
					local Ynew = e(ll)
					if `Ynew'<`target' {	// now bracketing target
						local left_limit = `bold'-`stepsize'*(`target'-`Yold')/(`Ynew'-`Yold')
					}
					else {
						local Yold `Ynew'
						local bold `A'
						local ++i
					}
					if "`dots'"!="nodots" noi di "." _c
				}
				local cost `i'
				if missing(`left_limit') noi di as txt _n "[note: failed to find left-hand confidence limit]"
			}
/*
		Now right limit
*/
			if `Y'[`n']>`target' {	// search for right_limit to the right of last value of b
				local Yold = `Y'[`n']
				local bold `to'
				local i 1
				while missing(`right_limit') & `i'<=`maxcost' {
					local A = `to'+`i'*`stepsize'
					// evaluate pll
					parsat `"`anything'"' if `touse', formula(`formula') var(`xx') value(`A') placeholder(`placeholder')
					cap `cmd' `r(result)' if `touse' `wt', `options'
					if _rc {
						noi di as err "model fit failed at parameter = " `A'
						exit 198
					}
					local Ynew = e(ll)
					if `Ynew'<`target' {	// now bracketing target
						local right_limit = `bold'+`stepsize'*(`target'-`Yold')/(`Ynew'-`Yold')
					}
					else {
						local Yold `Ynew'
						local bold `A'
						local ++i
					}
					if "`dots'"!="nodots" noi di "." _c
				}
				local cost = `cost'+`i'
			}
			else {	// pll is below target on right of mle, bracket target and interpolate
				local n1 = `n'-1
				forvalues i=`n1'(-1)1 {
					if `Y'[`i']>=`target' {
						local right_limit = `X'[`i']+`stepsize'*(`target'-`Y'[`i'])/(`Y'[`i'+1]-`Y'[`i'])
						continue, break	// exit forvalues loop
					}
				}
			}
			if missing(`right_limit') noi di as txt _n "[note: failed to find right-hand confidence limit]"
		}
		local se .
		local llci .
		local ulci .
		lab var `X' "`placeholder' in `formula'"
	}
}
// Pseudo-SE
local pse = (`right_limit'-`left_limit')/(2*`z')
capture drop `gen1'
capture drop `gen2'
rename `X' `gen1'
rename `Y' `gen2'
lab var `gen2' "profile log-likelihood function"
local ll_limit = `ll0'-`z'^2/2
local limit `ll_limit'
if "`difference'"!="" {
	// compute difference, subtract ll0
	qui replace `gen2' = `gen2'-`ll0'
	lab var `gen2' "profile log-likelihood difference function"
	local limit = -`z'^2/2
}
if "`deviance'"!="" {
	qui replace `gen2' = -2*`gen2'
	if "`difference'"!="" lab var `gen2' "profile deviance difference function"
	else lab var `gen2' "profile deviance function"
	local limit = -2*`limit'
}
local asym = 100*((`right_limit'-`b0')-(`b0'-`left_limit'))/(`right_limit'-`left_limit')
if "`graph'"!="nograph" {
	local Asym: display %4.1f `asym'

	// Extract title if present - default involves asymmetry
	_get_gropts , graphopts(`gropt') getallowed(title)
	local gropt `s(graphopts)'
	if `"`s(title)'"'=="" local title title("Asymmetry = `Asym'%")
	else if `"`s(title)'"'==`""""' local title
	else local title title(`s(title)')

	if !missing(`left_limit')  local lll `left_limit'
	if !missing(`right_limit') local rrr `right_limit'
	if ("`lll'`rrr'"!="")					///
		local xl xline(`lll' `rrr', lstyle(ci) `cilines')
	graph twoway line `gen2' `gen1', `gropt' `title' ///
	    `xl' yline(`limit', lstyle(refline) `levline')
}
if "`dots'"!="nodots" di
local tt "Coef."
di as txt _n "{hline 13}{c TT}{hline 47}"
di as txt %12s abbrev("`ytitle'",12) _col(14)"{c |}" ///
  %10s "Coef." "   Std. Err.     [`level'% PLL Conf. Int.]"
di as txt "{hline 13}{c +}{hline 47}"
if "`eq'"!="" di as res %-12s abbrev("`eq'",12) _col(14) as txt "{c |}"
if "`formula'"!="" di as txt %12s "`placeholder'" _cont
else di as txt %12s abbrev("`profile'",12) _cont
di _col(14) "{c |}" as res ///
  _col(16) %9.0g `b0' ///
  _col(28) %9.0g `pse'  ///
  _col(41) %9.0g `left_limit' ///
  _col(53) %9.0g `right_limit'
di as txt "{hline 13}{c BT}{hline 47}"
di as txt "Note: Std. Err. is pseudo standard error, derived from PLL CI"

// MLE of beta
return scalar b = `b0'
return scalar se = `se'
return scalar pse = `pse'

// Upper confidence limits
return scalar n_ulci = `ulci'
return scalar l_ulci = `right_limit'

// Lower confidence limits
return scalar n_llci = `llci'
return scalar l_llci = `left_limit'

// maximised likelihood
return scalar ll = `ll0'

// target ll for computing likelihood-based confidence limits
return scalar ll_limit = `ll_limit'
return scalar nobs = `nobs'

// Cost of searching for pll-based confidence interval
return scalar cost = `cost'

// Percentage asymmetry
return scalar asym = `asym'
end

program define parsat, rclass
version 9.0
syntax anything [if], Formula(string) var(varname) VALue(string) Placeholder(string)

// Replace placeholder with var and update var
local result: subinstr local anything "`placeholder'" " `var' ", all

// Substitute `placeholder' with `value' in `formula'
local f = subinstr(`"`formula'"', `"`placeholder'"' , "`value'", .)
quietly replace `var' = `f' `if'
return local result `result'
end
 
