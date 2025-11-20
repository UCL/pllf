/*
*! v1.3.7 PR 04mar2023 / IW 20nov2025
	better reporting of within-loop errors
	doesn't leave constraints lying around 
	advises replacing mixed with meglm 
v1.3.6 PR 04mar2023 / IW 04nov2025
	Syntax 2: if formula evaluates to constant, replace with its first derivative.
v1.3.5 PR 04mar2023 / IW 03nov2025
	drop collinearity check in syntax 2: correctly allows progress even when formula is constant
	new normcoll option
	corrected bug that made commands with offset() fail
	new warning if e(vce) is not "oim" or "eim"
v1.3.4 PR 04mar2023 / IW 12sep2025
	works for stcox with one covariate
	PLLF stored as double - matters for very large samples
	make trace work with syntax 2
v1.3.3 PR 04mar2023 / IW 02sep2025
	shownormal -> normal, and saved by gen()
	syntax 2: placeholder defaults to @ and is not needed on RHS
	don't ereturn results
	graph respects eform
v1.3.2 PR 04mar2023 / IW 01aug2025
	new shownormal option
	fix with stpm
	fix noconstant option
v1.3.1 PR 04mar2023 / IW 23jul2025
	mleline option
	new attempt to detect collinearity in xvars 
		fails by default, but new dropcollinear option forces continue
v1.3.0 PR 04mar2023 / IW 18jun2025
	sort out collinearity
	better error capturing
	allow profiling the constant
	new eform, trace and verbose options
v1.2.0 PR 04mar2023 / IW 28may2025	
	reworked as a prefix command
	fixed to work for prefect prediction

Code structure:
	parse pre-colon part
	parse post-colon part
	make checks across the parts
	run with profile() option (~200 lines)
	run with formula() option (~200 lines; calls parsat)
	compute pseudo SE, draw graph, print results, return results (~100 lines)
*/
program define pllf, rclass sortpreserve	
version 11.0

// SEPARATE PLLF/PREFIX PART FROM STATA COMMAND

mata: _parse_colon("hascolon", "statacmd")
if !`hascolon' {
	di as error "pllf is now a prefix command, with syntax:"
	di as error "    pllf<, options>: <regcmd>"
	exit 198
}

// PARSE PLLF OPTIONS FROM FIRST ARGUMENT

syntax [, ///
 DEViance FORMula(string) gen(string) DIFFerence LEVel(cilevel) ///
 PLaceholder(string) PROfile(string) range(string) ///
 MAXCost(int -1) N_eval(integer 100) noci noDOTs nograph gropt(string asis) ///
 LEVLINe(string asis) CILINes(string asis) ///
 TRace VERbose eform EFORM2(string) ///
 MLELine DROPCollinear NORMal NORMal2(string) /// to be documented
 debug debug2 LIst tol(real 1E-4) pause PROFILEOptions(string) noRMColl /// to remain undocumented
 ]

if `maxcost'<0 local maxcost = int(`n_eval'/2)

if "`placeholder'"!="" {
	if wordcount("`placeholder'")!=1 {
		di as err "invalid placeholder()"
		exit 198
	}
}
else local placeholder @

if "`gen'"!="" {
	local gen1 : word 1 of `gen'
	local gen2 : word 2 of `gen'
	local gen3 : word 3 of `gen'
}
if "`gen1'"=="" local gen1 _beta
if "`gen2'"=="" local gen2 _pll
if "`gen3'"=="" local gen3 _pllnorm

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

if !missing("`trace'") local dots nodots

if !missing("`debug'") local verbose verbose

if !missing("`verbose'") local noisily noisily
else local noisily quietly

if !missing("`normal2'") {
	local normal normal 
	local normalopts `normal2' 
}

*** END OF PARSING PLLF OPTIONS

*** PARSE REGRESSION COMMAND

gettoken cmd statacmd: statacmd

if substr("`cmd'", -1, .) == "," {
	local cmd = substr("`cmd'", 1, length("`cmd'") - 1)
	local statacmd ,`statacmd'
}

/*
	Currently supported commands are listed in help file and tested in test_regcmds.do
*/

else if ("`cmd'"=="fit") | ("`cmd'"=="reg") | (substr("`cmd'",1,4)=="regr") local cmd regress
if "`cmd'"=="reg3" {
	di as error "Sorry, reg3 is not supported"
	exit 498
}
if "`cmd'"=="mixed" {
	di as error "Sorry, mixed is not supported. Try meglm with the same syntax"
	exit 498
}

local 0 `statacmd'
syntax [anything] [if] [in] [using] [fweight pweight aweight iweight], [irr or hr coef NOHR offset(varname) exposure(varname) NOCONStant *]
if "`cmd'"=="logistic" & !missing("`coef'") local or or
if "`cmd'"=="stcox"  & !missing("`nohr'") local hr hr
if !missing("`irr'`or'`hr'") {
	local eform eform
	local eformname = upper("`irr'`or'`hr'")
}
if !missing("`eform2'") {
	local eform eform
	local eformname `eform2'
}
local constant `noconstant'
// Process user offset, if specified
if "`exposure'"!="" {
	if "`offset'"!="" {
		di as error "Can't specify both exposure() and offset()"
		exit 198
	}
	tempvar offset
	qui gen `offset' = log(`exposure')
}
if "`offset'"!="" {
	if "`cmd'"=="regress" {
		di as err "option offset() not allowed with regress"
		di as err "try first subtracting the offset `offset' from the outcome variable"
		exit 198
	}
	local plusoffset + `offset'
	local useroffset offset(`offset')
	local offset
	local niceuseroffset = cond(missing("`exposure'"), "`useroffset'", "exposure(`exposure')")
}

if "`weight'" != "" local wt [`weight'`exp']

if missing("`rmcoll'") & !missing("`anything'") { // -normcoll- not specified 
	unab varlist : `anything'
	* attempt to detect collinearity in xvars
	if inlist("`cmd'", "streg", "stcox", "stpm", "stpm2") {
		local yvar
		local xvars `varlist'
	}
	else gettoken yvar xvars : varlist
	if inlist("`cmd'","blogit") {
		gettoken yvar2 xvars : xvars
		local yvar `yvar' `yvar2'
	}
	_rmcoll `xvars' `if' `in' `wt', forcedrop `constant'
	if r(k_omitted)>0 {
		di as error "Collinearity found in xvarlist: `xvars'"
		if !missing("`dropcollinear'") {
			local xvars = r(varlist)
			di as error "dropcollinear option --> reduced xvarlist: `xvars'"
			local varlist `yvar' `xvars'
		}
		else if !missing("`allowcollinear'") di as text "proceeding with caution"
		else exit 498
	}
}
else local varlist `anything' // may include punctuation

*** END OF PARSING REGRESSION COMMAND

*** MIXED PARSING

if "`formula'"!="" {
	if "`profile'"!="" {
		di as txt "[profile() ignored]"
		local profile
	}
	if "`range'"=="" {
		di as err "range() required"
		exit 198
	}
}
else if "`profile'"!="" {
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
		if "`cmd'"=="stpm" {
			di as error "pllf cannot handle eqnames with stpm: please try stpm2"
			exit 498
		}
	}
	if "`profile'"!="_cons"{
		unab profile: `profile'
		local profilevar `profile'
	}
	else {
		local profilevar 1
	}
	* check `profile' is in `varlist'
	if "`profile'"=="_cons" & "`cmd'"=="stcox" {
		di as error "_cons does not exist in `cmd'"
		exit 198
	}
	if missing("`eq'") {
		if "`profile'"!="_cons" {
			local ok : list profile in varlist 
			if !`ok' {
				di as error "Profile variable `profile' not found in varlist"
				exit 198
			}
		}
	}
	if !missing("`debug'") {
		di as input "Debug: " _c
		if !missing("`eq'") di as input "equation is `eq'; " _c
		di as input "coefficient is `profile'"
	}
}
else {
	di as err "must supply either formula() or profile()"
	exit 198
}

*** END OF MIXED PARSING

*** START OF CODE FOR LINEAR PROFILING - PROFILE() OPTION
if "`profile'" != "" { // ------------ begin linear profiling --------
	// Fit model and get level% ci. Program terminates if invalid cmd attempted.
	if "`cmd'"=="stpm" local eq [xb]
	`noisily' di as text "Finding MLE using: " stritrim(`"`cmd' `varlist' `if' `in' `wt', `options' `constant' `niceuseroffset'"')
	capture `noisily' `cmd' `varlist' `if' `in' `wt', `options' `constant' `niceuseroffset'
	if _rc {
		di as error "Command failed: " stritrim(`"`cmd' `varlist' `if' `in' `wt', `options' `constant' `niceuseroffset'"')
		exit _rc
	}
	local ytitle `e(depvar)'
	local vcetype=e(vce)
	if !inlist("`vcetype'","oim","eim") di as text "Warning: pllf is a likelihood-based method and will ignore variance-covariance method `=e(vce)'"
	quietly {
		// Check that alleged parameter exists. 
		capture local b0 = `eq'_b[`profile']
		if _rc local b0 .
		capture local se = `eq'_se[`profile']
		if _rc local se .
		if `b0'==. di as error "Warning: after MLE, `eq'_b[`profile'] is missing"
		if `se'<=0 {
			di as error "Warning: after MLE, `eq'_se[`profile'] is zero" 
			di "`eq'_b[`profile'] and `eq'_se[`profile'] have been set to missing"
			local b0 .
		}
		if `se'==. {
			di as error "Warning: after MLE, `eq'_se[`profile'] is missing" 
			di "`eq'_b[`profile'] has also been set to missing"
		}
		if "`cmd'"=="stpm" local eq 

		if "`cmd'"=="regress" local z = invttail(e(df_r), (100-`level')/200)
		else local z = -invnorm((100-`level')/200)
		local nobs = e(N)
		local ll0  = e(ll)
		local use_deviance 0
		if missing(`ll0') {
			noi di as txt "note: valid log likelihood not returned in e(ll);"
			noi di as txt "if available, trying e(deviance) instead"
			local ll0 = -e(deviance)/2
			if missing(`ll0') {
				di as err "valid deviance not returned in e(deviance)"
				exit 198
			}
			else {
				local use_deviance 1
			}
		}
		if !missing("`se'") {
			local llci = `b0'-`z'*`se'
			local ulci = `b0'+`z'*`se'
		}
		else {
			if "`range'"=="" {
				noi di as err "Could not select range for `profile' (could not estimate MLE). Please supply range()"
				exit 198
			}
			local llci .
			local ulci .
		}
		if "`range'"=="" {
			// previously used default range as Wald-based `level' CI; now using +/-(z*1.2)SE for this.
			if `se'==. | `se'<=0 {
				di as error "Can't find SE, so can't choose a range. Please specify range()"
				exit 498
			}
			local from = `b0'-`z'*1.2*`se'
			local to =   `b0'+`z'*1.2*`se'
		}

		// If no equation specified, unabbreviate varlist and remove `profile' from it
		if "`eq'"=="" {
			// unab varlist: `varlist' // not needed, and fails with mixed
			if "`profile'"!="_cons" {
				local varlist: list varlist - profile
			}
			else {
				local constant noconstant
			}
		}

		if "`cmd'"=="regress" {
			constraint free
			local cuse `r(free)'

			// For regress, need to identify yvar; otherwise, not required
			gettoken yvar varlist: varlist
		}
		if `n_eval'>_N {
			di as txt "[Increasing the dataset size to `n_eval'] " _cont
			set obs `n_eval'
		}
		tempvar X Y order // values of X = regression coefficient, Y = profile likelihood for `profile'
		gen `X' = .
		gen double `Y' = .
		gen long `order' = _n
		local stepsize = (`to'-`from')/(`n_eval'-1)
		if !missing("`debug'") noisily di
		noisily di as text "Profiling:" _c
		local anyerror1 0
		local anyerror2 0
		forvalues i=1/`n_eval' { // MAIN LOOP FOR PROFILE()
			ereturn clear // to avoid picking up old results
			local b = `from'+(`i'-1)*`stepsize'
			if "`eq'"!="" {
				// Use constrained regression
				if !missing("`debug'") & `i'==1 noi di as text " using constraint" _c
				constraint define `cuse' `eq'`profile'=`b'
				local thiscmd `cmd' `varlist' `if' `in' `wt', `options' constraint(`cuse') `constant' `useroffset' `profileoptions'
			}
			else {
				if `i'==1 {
					tempvar offset
					gen `offset' = .
				}
				if "`cmd'"=="regress" {
					if !missing("`debug'") & `i'==1 noi di as text " using modified outcome" _c
					replace `offset' = `yvar' - `b'*`profilevar'
					local thiscmd regress `offset' `varlist' `if' `in' `wt', `options' `constant' `profileoptions'
				}
				else {
					if !missing("`debug'") & `i'==1 noi di as text " using offset" _c
					replace `offset' = `b'*`profilevar' `plusoffset'
					if "`cmd'"=="stcox" & missing("`varlist'") local estimate estimate
					local thiscmd `cmd' `varlist' `if' `in' `wt', `options' offset(`offset') `constant' `estimate' `profileoptions'
				}
			}
			if mi("`debug'") cap `thiscmd'
			else {
				noi di as input _new(2) "Parm: `profile' = `b'"
				noi di `"Command is: `thiscmd'"'
				if !mi("`cuse'") noi constraint dir `cuse'
				cap noi `thiscmd'
			}
			if _rc==1 exit 1
			if _rc==198 {
				noi di as error "pllf error: syntax error when fixing parameter value using:"
				noi di as error `"`thiscmd'"'				
				exit 198
			}
			local error = _rc
			if _rc==0 & e(converged)==0 local error 1
			if !missing("`debug2'") {
				noisily di as input "coeff=`b'"
				noisily `cmd'
			}
			if !missing("`trace'") & `i'==1 noi di ":"
			sort `order'
			replace `X' = `b' in `i'
			replace `Y' = cond(`use_deviance', -e(deviance)/2, e(ll)) in `i'
			if "`dots'"!="nodots" {
				if `error'>1 {
					noi di as error "x[`error']" _c
					local anyerror2 1
				}
				else if `error'==1 {
					noi di as error "?" _c
					local anyerror1 1
				}
				else noi di as text "." _c
			}
			if !missing("`trace'") {
				local col = length("`profile'")+14
				noi di as text "`profile'=" as result `b' _c _col(`col')
				local resultfmt = cond(`error',"error","result")
				if !`use_deviance' noi di as text "PLL=" as `resultfmt' e(ll) _c
				else noi di as text "deviance=" as `resultfmt' e(deviance) _c
				if `error'>2 noi di " [error `error']"
				noi di
			}
		}
		if `anyerror1' noi di as error _n "?" as text " means `cmd' did not converge for this parameter value"
		if `anyerror2' noi di as error _n "x" as text " means `cmd' failed for this parameter value"
		if !`anyerror1' & !`anyerror2' noi di
		summ `Y', meanonly
		if r(max)-r(min) < `tol' {
			di as error "Warning: PLL does not seem to vary over this range"
			exit 498
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
				forvalues i=2/`n_eval' {
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
				noisily di as text "Lower confidence limit not yet found: searching" _c
				if !missing("`trace'") noi di
				local Yold = `Y'[1]
				local bold `from'
				local i 1
				while missing(`left_limit') & `i'<=`maxcost' {
					local b = `from'-`i'*`stepsize'
					// evaluate pll
					if "`eq'"!="" {
						// Use constrained regression
						constraint define `cuse' `eq'`profile'=`b'
						`cmd' `varlist' `if' `in' `wt', `options' constraint(`cuse') `constant' `useroffset' `profileoptions'
					}
					else {
						if "`cmd'"=="regress" {
							replace `offset' = `yvar' - `b'*`profilevar'
							regress `offset' `varlist' `if' `in' `wt', `options'
						}
						else {
							replace `offset' = `b'*`profile' `plusoffset'
							`cmd' `varlist' `if' `in' `wt', `options' offset(`offset') `constant'
						}
					}
					local Ynew = cond(`use_deviance', -e(deviance)/2, e(ll))
					if `Ynew'<`target' {	// now bracketing target
						local left_limit = `bold'-`stepsize'*(`target'-`Yold')/(`Ynew'-`Yold')
					}
					else {
						local Yold `Ynew'
						local bold `b'
						local ++i
					}
					if "`dots'"!="nodots" noi di as text "." _c
					else if !missing("`trace'") {
						local col = length("`profile'")+14
						noi di as text "`profile'=" as result `b' _c _col(`col')
						if !`use_deviance' noi di as text ", PLL=" as result `Ynew'
						else noi di as text ", deviance=" as result `Ynew'
					}

				}
				local cost `i'
				if missing("`trace'") noi di
				if missing(`left_limit') noi di as txt "[note: failed to find lower confidence limit]"
			}
/*
		Now right limit
*/
			if `Y'[`n_eval']>`target' {	// search for right_limit to the right of last value of b
				noisily di as text "Upper confidence limit not yet found: searching" _c
				if !missing("`trace'") noi di
				local Yold = `Y'[`n_eval']
				local bold `to'
				local i 1
				while missing(`right_limit') & `i'<=`maxcost' {
					local b = `to'+`i'*`stepsize'
					// evaluate pll
					if "`eq'"!="" {
						// Use constrained regression
						constraint define `cuse' `eq'`profile'=`b'
						`cmd' `varlist' `if' `in' `wt', `options' constraint(`cuse') `constant' `useroffset'
					}
					else {
						if "`cmd'"=="regress" {
							replace `offset' = `yvar' - `b'*`profile'
							regress `offset' `varlist' `if' `in' `wt', `options' `constant'
						}
						else {
							replace `offset' = `b'*`profile' `plusoffset'
							`cmd' `varlist' `if' `in' `wt', `options' offset(`offset') `constant'
						}
					}
					local Ynew = cond(`use_deviance', -e(deviance)/2, e(ll))
					if `Ynew'<`target' {	// now bracketing target
						local right_limit = `bold'+`stepsize'*(`target'-`Yold')/(`Ynew'-`Yold')
					}
					else {
						local Yold `Ynew'
						local bold `b'
						local ++i
					}
					if "`dots'"!="nodots" noi di as text "." _c
					else if !missing("`trace'") {
						noi di as text "`profile'=" as result `b' _c
						if !`use_deviance' noi di as text ", PLL=" as result `Ynew'
						else noi di as text ", deviance=" as result `Ynew'
					}
				}
				local cost = `cost'+`i'
			}
			else {	// pll is below target on right of mle, bracket target and interpolate
				local n1 = `n_eval'-1
				forvalues i=`n1'(-1)1 {
					if `Y'[`i']>=`target' {
						local right_limit = `X'[`i']+`stepsize'*(`target'-`Y'[`i'])/(`Y'[`i'+1]-`Y'[`i'])
						continue, break	// exit forvalues loop
					}
				}
			}
			if missing("`trace'") noi di
			if missing(`right_limit') noi di as txt "[note: failed to find upper confidence limit]"
		}
		if missing("`eform'") lab var `X' "`eq'_b[`profile']"
		else lab var `X' "exp(`eq'_b[`profile'])"
	}
}
cap constraint drop `cuse'

*** END OF CODE FOR LINEAR PROFILING
*** START OF CODE FOR NON-LINEAR PROFILING - FORMULA() OPTION
else {	// --------------- begin non-linear profiling ---------------
	tempvar xx
	qui gen `xx' = .
	// Trial fit of model with central value of param. Program terminates if invalid cmd attempted.
	local A = (`to'+`from')/2 // IW 
	local stepsize = (`to'-`from')/(`n_eval'-1)
	local varlist `varlist' `placeholder'
	parsat `"`varlist'"' `if' `in', formula(`formula') var(`xx') value(`A') placeholder(`placeholder') step(`stepsize')
	local result `r(result)'
	CheckCollin `result'
	if ("`s(result)'" == "") {
		local A = `A' * 1.05 // adjust a bit
		parsat `"`varlist'"' `if' `in', formula(`formula') var(`xx') value(`A') placeholder(`placeholder') step(`stepsize')
	}
	if !missing("`debug'") di as input "Starting with parm = `A'"
	if !missing("`debug'") di as input `"MLE: `cmd' `result' `if' `in' `wt', `options' `constant' `useroffset'"'
	cap `cmd' `result' `if' `in' `wt', `options' `constant' `useroffset'
	local rc = _rc
	if `rc' error `rc'

	local ytitle `e(depvar)'
	local ll = e(ll)
	quietly {
		if "`cmd'"=="regress" local z = invttail(e(df_r), (100-`level')/200)
		else local z = -invnorm((100-`level')/200)

		local nobs = e(N)

		if `n_eval'>_N {
			di as txt "[Increasing the dataset size to `n_eval'] " _cont
			set obs `n_eval'
		}
		tempvar X Y order // values of X = regression coefficient, Y = profile likelihood for `profile'
		gen `X' = .
		gen double `Y' = .
		gen long `order' = _n
		local y1 .
		local y2 .
		local y3 .
		local b0 .
		local ll0 .
		local done 0
		forvalues i=1/`n_eval' { // MAIN LOOP FOR FORMULA()
			local A = `from'+(`i'-1)*`stepsize'
			// Substitute A in `formula' and fit model
			parsat `"`varlist'"' `if' `in', formula(`formula') var(`xx') value(`A') placeholder(`placeholder') step(`stepsize')
			// Search for collinearity in expression, e.g. caused by x and x^p with p = 1.
			local result `r(result)'
summ `xx'
if r(sd)==0 noi di as error "formula gives SD=0 for parm=`A'"
			CheckCollin `result'
			if "`s(result)'" == "" {
				noi di as err "collinearity detected with parameter value = " `A'
				noi di as err "you may need to exclude this value from the parameter range"
				*exit 198
			}
			cap `noisily' `cmd' `result' `if' `in' `wt', `options' `constant' `useroffset'
			if _rc {
				noi di as err "model fit failed at parameter = " `A'
				exit 198
			}
			local vcetype=e(vce)
			if `i'==1 & !inlist("`vcetype'","oim","eim") di as error "{p 0 2}Warning: pllf is a likelihood-based method and will ignore the variance-covariance method calculated by method `=e(vce)'{p_end}"
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
			if "`dots'"!="nodots" noi di as text "." _c
			if !missing("`trace'") {
				local col 25
				noi di as text "Parameter `placeholder'=" as result `A' _c _col(`col')
				noi di as text "PLL=" as result e(ll)
			}
		}
		if !`done' {
			// MLE not straddled. Attempt quadratic solution using terminals of range and a midpoint. Maintain equal spacing.
			if mod(`n_eval',2)==0 {	// even number of evaluations
				local mid = `n_eval'/2
				local last `n_eval'-1
			}
			else {	// odd number of evaluations
				local mid = (`n_eval'+1)/2
				local last `n_eval'
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
				forvalues i=2/`n_eval' {
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
					parsat `"`varlist'"' `if' `in', formula(`formula') var(`xx') value(`A') placeholder(`placeholder') step(`stepsize')
					cap `cmd' `r(result)' `if' `in' `wt', `options' `constant' `useroffset'
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
					if "`dots'"!="nodots" noi di as text "." _c
				}
				local cost `i'
				if missing(`left_limit') noi di as txt _n "[note: failed to find left-hand confidence limit]"
			}
/*
		Now right limit
*/
			if `Y'[`n_eval']>`target' {	// search for right_limit to the right of last value of b
				local Yold = `Y'[`n_eval']
				local bold `to'
				local i 1
				while missing(`right_limit') & `i'<=`maxcost' {
					local A = `to'+`i'*`stepsize'
					// evaluate pll
					parsat `"`varlist'"' `if' `in', formula(`formula') var(`xx') value(`A') placeholder(`placeholder') step(`stepsize')
					cap `cmd' `r(result)' `if' `in' `wt', `options' `constant' `useroffset'
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
					if "`dots'"!="nodots" noi di as text "." _c
				}
				local cost = `cost'+`i'
			}
			else {	// pll is below target on right of mle, bracket target and interpolate
				local n1 = `n_eval'-1
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
		if missing("`eform'") lab var `X' "`placeholder' in `formula'"
		else lab var `X' "exp(`placeholder') in `formula'"
		if "`dots'"!="nodots" noi di
	}
	local use_deviance 0
}
*** END OF CODE FOR NON-LINEAR PROFILING

*** FINAL CODE COMMON TO BOTH PROFILE AND FORMULA

if !missing("`list'") l `X' `Y' if !missing(`X')

// Pseudo-SE
local pse = (`right_limit'-`left_limit')/(2*`z')

if !missing("`normal'") {
	if !missing("`profile'") local usese `se'
	else local usese `pse'
	if missing(`usese') | `usese'<=0 local usese `pse'
	if missing(`usese') | `usese'<=0 {
		noisily di as error "Warning: can't graph Normal approximation when se is missing"
		local normal
	}
	else {
		tempvar Z
		qui gen double `Z' = -0.5*((`X'-`b0')/`usese')^2
		lab var `Z' "normal approximation"
	}
}

capture drop `gen1'
rename `X' `gen1'
capture drop `gen2'
rename `Y' `gen2'
if !missing("`normal'") {
	capture drop `gen3'
	rename `Z' `gen3'
}
if `use_deviance' local star = "*"
lab var `gen2' "profile log likelihood`star'"
local ll_limit = `ll0'-`z'^2/2
local limit `ll_limit'
if "`difference'"!="" {
	// compute difference, subtract ll0
	qui replace `gen2' = `gen2'-`ll0'
	lab var `gen2' "profile log likelihood difference`star'"
	local limit = -`z'^2/2
}
else if !mi("`normal'") qui replace `gen3' = `gen3'+`ll0'
if "`deviance'"!="" {
	qui replace `gen2' = -2*`gen2'
	if !mi("`normal'") qui replace `gen3' = -2*`gen3'
	if "`difference'"!="" lab var `gen2' "profile deviance difference`star'"
	else lab var `gen2' "profile deviance`star'"
	local limit = -2*`limit'
}
if !missing("`mleline'") local limit `limit' `ll0'
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
	if !missing("`eform'") {
		qui replace `gen1'=exp(`gen1') // will be reversed later
		if !missing("`lll'") local lll = exp(`lll')
		if !missing("`rrr'") local rrr = exp(`rrr')
		cap niceloglabels `gen1', local(betalabel) style(125)
		if !_rc local gropt xscale(log) xlabel(`betalabel') `gropt'
		if _rc==199 di as text as smcl "To improve xlabels, please use {stata ssc install niceloglabels}"
	}
	if ("`lll'`rrr'"!="")					///
		local xl xline(`lll' `rrr', lstyle(ci) `cilines')
	if `use_deviance' {
		local note note("`star'Defining log likelihood = -0.5*e(deviance)")
	}
	if !mi("`normal'") {
		local normgraph (line `gen3' `gen1', lpattern(dash) `normalopts')
		local ytitle : variable label `gen2'
		local gropt ytitle(`ytitle') `gropt'
	}
	local graphcmd graph twoway (line `gen2' `gen1') `normgraph', `gropt' `title' ///
	    `xl' yline(`limit', lstyle(refline) `levline')
	if !missing("`debug'") di as input `"Drawing graph: `graphcmd'"'
	if !missing("`pause'") {
		global F9 `graphcmd'
		pause
	}
	`graphcmd'
	if !missing("`eform'") qui replace `gen1'=log(`gen1')
}
di as txt _n "{hline 13}{c TT}{hline 47}"
di as txt %12s abbrev("`ytitle'",12) _col(14)"{c |}" _c
local parmname = cond(!missing("`profile'"),"Coef","Est")
if missing("`eform'") di as txt %10s "`parmname'." _c
else if !missing("`eformname'") di as txt %10s "`eformname'" _c
else di as txt %10s "exp(parmname)" _c
di as txt "   Std. Err.     [`level'% PLL Conf. Int.]"
di as txt "{hline 13}{c +}{hline 47}"
if "`eq'"!="" di as res %-12s abbrev("`eq'",12) _col(14) as txt "{c |}"
if "`formula'"!="" di as txt %12s "Parameter `placeholder'" _cont
else di as txt %12s abbrev("`profile'",12) _cont
if missing("`eform'") {
	local b0_display `b0'
	local pse_display `pse'
	local left_limit_display `left_limit'
	local right_limit_display `right_limit'
}
else {
	local b0_display exp(`b0')
	local pse_display `pse'*exp(`b0')
	local left_limit_display exp(`left_limit')
	local right_limit_display exp(`right_limit')
}
di _col(14) "{c |}" as res ///
  _col(16) %9.0g `b0_display' ///
  _col(28) %9.0g `pse_display'  ///
  _col(41) %9.0g `left_limit_display' ///
  _col(53) %9.0g `right_limit_display'
di as txt "{hline 13}{c BT}{hline 47}"
di as txt "Note: Std. Err. is pseudo standard error, derived from PLL CI"

*** RETURN RESULTS

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

// Don't keep last model in memory, since it isn't the MLE
ereturn clear

end




program define parsat, rclass
version 9.0
syntax anything [if], Formula(string) var(varname) VALue(string) Placeholder(string) step(real)

// Replace placeholder with var and update var
local result: subinstr local anything "`placeholder'" " `var' ", all

// Substitute `placeholder' with `value' in `formula'
local f = subinstr(`"`formula'"', `"`placeholder'"' , "`value'", .)
quietly replace `var' = `f' `if'

// check if no variance
quietly summ `var'
if r(sd)==0 { // replace with first (numerical) derivative
	local eps = `step'/1000
	local vl = `value'-`eps'
	local fl = subinstr(`"`formula'"', `"`placeholder'"' , "`vl'", .)
	local vu = `value'+`eps'
	local fu = subinstr(`"`formula'"', `"`placeholder'"' , "`vu'", .)
	quietly replace `var' = ((`fu')-(`fl'))/(2*`eps') `if'	
}

return local result `result'
end




program define CheckCollin, sclass
version 9.0
local l1 : word count `*'
quietly _rmcoll `*'
local result `r(varlist)'
local l2 : word count `result'
if (`l2' < `l1') sreturn local result ""
else sreturn local result `result'
end
