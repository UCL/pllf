/*
Run all tests for pllf
IW 3nov2025
*/

* personal setup
global path c:\ian\git\pllf
* end of personal setup

adopath ++ $path
prog drop _all
cap log close
foreach filename in test_pllf test_regcmds test_new_rmcoll TRISST {
	cd $path\testing
	set linesize 100
	clear all 

	log using `filename', replace text
	version
	which pllf
	cap noi do `filename'
	log close
	if _rc local failfile `failfile' `filename'
}

if mi("`failfile'") di "pllf passed all its tests"
else di as error "pllf failed in: `failfile'"
