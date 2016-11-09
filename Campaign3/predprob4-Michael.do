// A workflow for two- or three-level logistic regression models with cross-level interactions
// Frederick Solt, October 2012

/*************************************************
This workflow:

(1) generates estimates for two- or three-level  
logistic regressions of any number of dependent  
variables (but ideally three or fewer) that include
interactions (i.e., random effects) between any 
number of specified level-1 variables, on the one hand, 
and specified level-2 variables, on the other hand; 
the model must be identical for all DVs, and all 
level-1 interactions must be with the same set of
level-2 variables;

(2) generates predicted probabilities from the 
estimates, which are captured in separate log files
for each DV.

(3) generates graphs of the predicted probabilities
over the range of each level-2 variable for selected
values (for continuous variables, the minimum, median,
and maximum) of each level-1 variable.

(4) generates a dotplot of the logit results, with
coefficients of continuous/ordinal IVs rescaled by 
multiplying by 2 s.d.s so as to be on a similar scale
to the (unchanged) coefficients of dichotomous IVs.
See Kastellec and Leoni (2007) regarding dotplots, 
and Gelman (2008) regarding rescaling by 2 s.d.s.

Although this workflow _could_ be wrapped into a 
command fairly easily, this format invites users 
to make modifications to suit their needs.

Dependencies:
HLM 6.0 (and, on OS X, Crossover)
FS updates to Sean Reardon's hlm commands
distinct (by Gary Longton and Nicholas Cox, automatically installed if missing)
eclplot and sencode (by Roger Newson, automatically installed if missing)
**************************************************/

set more off
clear
capture log close
cd "~/Documents/Projects/Workflow/Michael/Campaign3"	/*Replace with path to your folder for analyses for this project*/
log using wf.log, replace				/*Replace wf with the name of your project*/

// First, describe data and model
local dataset = "cces3"			/*Replace data_set with name of your dataset that includes all variables*/

local dvs = "donate meetings polwork" 					/*Dependent variables*/
local l3_id = ""				/*Level-3 identification variable*/
local l2_id = "State2"					/*Level-2 identification variable*/
local iv_l3 = ""		/*Level-3 predictors*/
local iv_l2 = "gini_st1 gdppc div_st inits_st gov_st sen_st sw_mar" 		/*Level-2 predictors*/
local iv_l2_interact = "gini_st1"	/*Level-2 predictor(s) that interact with level-1 predictors*/
local iv_l1 = "age educated rep dem black hisp male married children church unionmem income" /*Level-1 predictors without interactions*/
local iv_l1_interact = "income"			/*Level-1 predictor(s) with level-2 interactions*/

// No further user modification _necessary_

// Confirm all needed Stata commands installed
capture ssc install distinct
capture ssc install eclplot
capture ssc install sencode

// Proceed
use `dataset'.dta, clear

if "`l3_id'"=="" & "`iv_l3'"=="" local lev = 2
else if "`l3_id'"=="" {
	di in red "When including level-3 predictors, the level-3 identification variable must be specified."
	error
}
else local lev = 3

local iv_l1_no_interact = regexr("`iv_l1'", "`iv_l1_interact'", "")
local iv_l1_int_model = ""
local iv_l1_all = "`iv_l1_no_interact'"

foreach a in `iv_l1_interact' {
	local iv_l1_add = "`a'"
	foreach b in `iv_l2_interact' {
		local iv_l1_add = "`iv_l1_add'" + " `a'X`b'"
	}
	local iv_l1_all = "`iv_l1_all'" + " `iv_l1_add'"
	local iv_l1_int_model = "`iv_l1_int_model'"+" `a'(int `iv_l2_interact' rand)"
}

local ivs = "cons `iv_l3' `iv_l2' `iv_l1_all'"
local total_ivs = wordcount("`ivs'")		
local total_dvs = wordcount("`dvs'")
forvalues i = 1/`total_ivs' {
	local t = word("`ivs'", `i')
	local `t'_place = `i'
}

local n_l1_no_i = wordcount("`iv_l1_no_interact'")
local n_l2_no_i = wordcount(regexr("`iv_l2'","`iv_l2_interact'", ""))
local n_l3 = wordcount("`iv_l3'")

gen cons=1
foreach a in `iv_l1_interact' {
	foreach b in `iv_l2_interact' {
	gen `a'X`b' = `a'*`b'
	}
}
foreach v in `dvs' `ivs'{
		local lab_`v' : variable label `v'
		if regexm("`v'","([^X]*)X([^X]*)") {
			local aa = regexs(1)
			local bb = regexs(2)
			if "`lab_`aa''" == "Income Quintile" & "`lab_`bb''" == "Income Inequality" {
				local lab_`v' = "Inequality*Income"
			}
			else if "`lab_`aa''" == "Income Quintile" {
				local lab_`v' = "`lab_`bb''"+"*"+"Income"
			}
			else if "`lab_`bb''" == "Income Inequality" {
				local lab_`v' = "Inequality"+"*"+"`lab_`aa''"
			}
			else local lab_`v' = "`lab_`bb''"+"*"+"`lab_`aa''"	
		}
		if "`lab_`v''"=="" local lab_`v' = proper("`v'")
}

local stub "b_"
tempfile ds1
save `ds1'

foreach dv in `dvs' {
	use `ds1'
	quietly levelsof `dv'
	if "`r(levels)'"=="0 1" {
		local nc = 2
		local model = "bernouli"
	}
	else {
		quietly distinct `dv'
		local nc = `r(ndistinct)'
		if `nc' == 2 {
			levelsof `dv'
			local val0 = word("`r(levels)'", 1)
			local val1 = word("`r(levels)'", 2)
			gen `dv'_2 = `dv'
			recode `dv'_2 `val0' = 0 `val1' = 1
			local dv = `dv'_2
		}
		else if `nc' <= 5 {
			local model = "ordinal ncat(`nc')"
			local tholds = "th2"
			if `nc' == 4 local tholds = "th2 th3"
			if `nc' == 5 local tholds = "th2 th3 th4"
			local total_tholds = wordcount("`tholds'")
			forvalues i = 1/`total_tholds' {
				local t = word("`tholds'", `i')
				local `t'_place = `total_ivs'+`i'
			}
		}	
		else {
			di in red "Dependent variable has too many categories to be treated as ordinal by HLM" 
			error
		}
	}
	if `nc' == 2 local sto = `"store("pa")"'
	
	use `dataset'.dta, clear	
	capture mkdir `dv'
	quietly cd `dv'
	
	if `lev'==2 {
		hlmmkmdm2 using `dv', id2(`l2_id') dv(`dv') l1(`iv_l1_no_interact' `iv_l1_interact') ///
			l2(`iv_l2') replace
		hlm hlm2 `dv' int(int `iv_l2' rand)  `iv_l1_no_interact' `iv_l1_int_model' rand, cmd("`dv'") ///
			mdmfile("`dv'.mdm") `model' `sto' robust out("`dv'.txt") replace run stop
	}
	else {
		hlmmkmdm3 using `dv', id2(`l2_id') id3(`l3_id') dv(`dv') l1(`iv_l1_no_interact' `iv_l1_interact') ///
			l2(`iv_l2') l3(`iv_l3') replace
		hlm hlm3 `dv' int(int(int `iv_l3' rand) `iv_l2' rand)  `iv_l1_no_interact' `iv_l1_int_model' rand, cmd("`dv'") ///
			mdmfile("`dv'.mdm") `model' `sto' robust out("`dv'.txt") replace run stop
	}
	use "`dv'_l1.dta", clear
	sort `l2_id'
	merge `l2_id' using "`dv'_l2.dta", nokeep
	drop _merge
	if `lev'==3 {
		sort `l3_id' `l2_id'
		merge `l3_id' using "`dv'_l3.dta", nokeep
		drop _merge
	}
	gen cons=1
	if `nc' > 2 {
		foreach i in `tholds' {
			gen `i' = 1
		}
	}
	foreach a in `iv_l1_interact' {
		foreach b in `iv_l2_interact' {
		gen `a'X`b' = `a'*`b'
		}
	}	

	set seed 324
	estsimp2, sims(2500) genname(`stub')
	save "`dv'_1.dta", replace
	
	//Create variables to store mean values of all variables	
	local i=0

	*Level-3 vars
	if `lev'==3 {
		use "`dv'_l3.dta", clear
		gen cons=1
		
		foreach a of varlist cons `iv_l3' {
			local ++i
			sum `a'
			gen v_`i'=r(mean)
		}
		
		sort `l3_id'
		save "`dv'_3.dta", replace	
	}
	else local cc = "cons"
	
	*Level-2 vars
	use "`dv'_l2.dta", clear
	
	if `lev'==2 gen cons=1
	
	foreach a of varlist `cc' `iv_l2' {
		local ++i
		sum `a'
		gen v_`i'=r(mean)
	}
	
	sort `l2_id'
	save "`dv'_2.dta", replace
	
	*Level-1 vars
	use "`dv'_1.dta", clear

	if `lev'==3 {
		sort `l3_id'
		merge `l3_id' using "`dv'_3.dta", _merge(_m3)
	}
	
	sort `l2_id'
	merge `l2_id' using "`dv'_2.dta", _merge(_m2)
	drop _m*
	
	foreach a of varlist `iv_l1_all' `tholds' {
		local ++i
		sum `a'
		gen v_`i'=r(mean)
	}
	
	save "`dv'_res.dta", replace

	//Generate baseline
	gen base=0
	forvalues i=1/`total_ivs' {
		replace base=base+(v_`i' * `stub'`i')
	}
	
	save "`dv'_res.dta", replace
	erase "`dv'_2.dta" 
	capture erase "`dv'_3.dta"
		
	//table
	capture log close `dv'_results
	log using "`dv'_results.log", replace name("`dv'_results")

	foreach a of varlist `iv_l1_interact' {
		foreach b of varlist `iv_l2_interact' {
			qui tab `a'
			if r(r) <= 5 {
				qui levelsof `a', local(`a'_vals)
			}
			else {
				qui sum `a'
				local min_`a' = round(r(min), .01)
				local max_`a' = round(r(max), .01)
				qui centile `a'
				local med_`a' = round(r(c_1), .01)
				local `a'_vals = "`min_`a'' `med_`a'' `max_`a''"
			}

			qui sum `b'
			local min_`b' = r(min)
			local max_`b' = r(max)
			
			local c = regexr("`iv_l2_interact'","`b'","")
			local c = trim("`c'")
			local total_others = wordcount("`c'")
			
			local ii = 1
			foreach i2 in ``a'_vals' {
				capture drop est_others
				gen sub_others = 0
				gen add_others = 0
				forvalues oth = 1/`total_others' {
					local other = word("`c'",`oth')
					qui replace sub_others = sub_others + v_``other'_place'*(`stub'``other'_place' ///
						+ v_``a'_place'*`stub'``a'X`other'_place')					
					qui replace add_others = add_others + v_``other'_place'*(`stub'``other'_place' ///
						+ `i2'*`stub'``a'X`other'_place')
				}
				qui gen est1`b'lo_`a'`ii' = base  ///
					- (v_``b'_place'*`stub'``b'_place'+ v_``a'_place'*`stub'``a'_place' ///
						+v_``a'_place'*v_``b'_place'*`stub'``a'X`b'_place'+sub_others)  ///
					+ (`min_`b''*`stub'``b'_place'+`i2'*`stub'``a'_place' ///
						+`i2'*`min_`b''*`stub'``a'X`b'_place'+add_others)
				qui gen est1`b'hi_`a'`ii' = base  ///
					- (v_``b'_place'*`stub'``b'_place'+ v_``a'_place'*`stub'``a'_place' ///
						+v_``a'_place'*v_``b'_place'*`stub'``a'X`b'_place'+sub_others)  ///
					+ (`max_`b''*`stub'``b'_place'+`i2'*`stub'``a'_place' ///
						+`i2'*`max_`b''*`stub'``a'X`b'_place'+add_others)
						
				qui gen pp1`b'lo_`a'`ii'=1/(1+exp(-est1`b'lo_`a'`ii'))
				label var pp1`b'lo_`a'`ii' "Predicted prob of `dv'=1 at min `b' when `a' = `i2'"
				qui gen pp1`b'hi_`a'`ii'=1/(1+exp(-est1`b'hi_`a'`ii')) 
				label var pp1`b'hi_`a'`ii' "Predicted prob of `dv'=1 at max `b' when `a' = `i2'"

				qui gen diff1`b'_`a'`ii'=pp1`b'hi_`a'`ii'-pp1`b'lo_`a'`ii'
				label var diff1`b'_`a'`ii' "Change in predicted prob of `dv'=1 over range of `b' when `a' = `i2'"

				di in white _newline(2) "for `a' = " `i2'
				sum pp1`b'lo_`a'`ii' pp1`b'hi_`a'`ii' diff1`b'_`a'`ii'
				
				if `nc' > 2 {
					qui gen tot_prob_lo = pp1`b'lo_`a'`ii'
					qui gen tot_prob_hi = pp1`b'hi_`a'`ii'
	
					forvalues c = 2/`nc' {
						if `c' == `nc' {
							qui gen pp`c'`b'lo_`a'`ii' = 1-tot_prob_lo
							label var pp`c'`b'lo_`a'`ii' "Predicted prob of `dv'=`c' at min `b' when `a' = `i2'"

							qui gen pp`c'`b'hi_`a'`ii' = 1-tot_prob_hi
							label var pp`c'`b'hi_`a'`ii' "Predicted prob of `dv'=`c' at max `b' when `a' = `i2'"							

							qui gen diff`c'`b'_`a'`ii'=pp`c'`b'hi_`a'`ii'-pp`c'`b'lo_`a'`ii'
							label var diff`c'`b'_`a'`ii' "Change in predicted prob of `dv'=`c' over range of `b' when `a' = `i2'"
						}
						else {
							qui gen est`c'`b'lo_`a'`ii' = est1`b'lo_`a'`ii'+(v_`th`c'_place'*`stub'`th`c'_place')
							qui gen pp`c'`b'lo_`a'`ii'=1/(1+exp(-est`c'`b'lo_`a'`ii'))-tot_prob_lo
							label var pp`c'`b'lo_`a'`ii' "Predicted prob of `dv'=`c' at min `b' when `a' = `i2'"
							replace tot_prob_lo = tot_prob_lo + pp`c'`b'lo_`a'`ii'
							
							qui gen est`c'`b'hi_`a'`ii' = est1`b'hi_`a'`ii'+(v_`th`c'_place'*`stub'`th`c'_place')
							qui gen pp`c'`b'hi_`a'`ii'=1/(1+exp(-est`c'`b'hi_`a'`ii'))-tot_prob_hi
							label var pp`c'`b'hi_`a'`ii' "Predicted prob of `dv'=`c' at max `b' when `a' = `i2'"							
							replace tot_prob_hi = tot_prob_hi + pp`c'`b'hi_`a'`ii'
							
							qui gen diff`c'`b'_`a'`ii'=pp`c'`b'hi_`a'`ii'-pp`c'`b'lo_`a'`ii'
							label var diff`c'`b'_`a'`ii' "Change in predicted prob of `dv'=`c' over range of `b' when `a' = `i2'"
						}
						sum pp`c'`b'lo_`a'`ii' pp`c'`b'hi_`a'`ii' diff`c'`b'_`a'`ii'
					}
				}			
				drop sub_others add_others
				if `nc' > 2 drop tot_prob_lo tot_prob_hi
				local ++ii
			}
		}
	}
			
	// Pred probs of other IVs, for comparison
	local j=0
	foreach iv in `ivs' {
		local j=`j'+1
		if regexm("`iv_l3' `iv_l2' `iv_l1_no_interact'", "`iv'")  ///
			& !(regexm("`iv_l2_interact'", "`iv'")) {
			qui sum `iv'
			local min = r(min)
			local max = r(max)
			
			qui gen est1lo_`j'=base-v_`j'*`stub'`j'+`min'*`stub'`j'
			qui gen est1hi_`j'=base-v_`j'*`stub'`j'+`max'*`stub'`j'			
			
			qui gen pp1lo_`j'=1/(1+exp(-est1lo_`j'))
			label var pp1lo_`j' "Predicted probability of `dv'=1 at min `iv'"			
			qui gen pp1hi_`j'=1/(1+exp(-est1hi_`j'))
			label var pp1hi_`j' "Predicted probability of `dv'=1 at max `iv'"
	
			qui gen diff1pp_`j'=pp1hi_`j'-pp1lo_`j'
			label var diff1pp_`j' "Change in predicted probability of `dv' over range of `iv'"
			
			di in white "`iv'" ", " `min' " to " `max'
			sum pp1lo_`j' pp1hi_`j' diff1pp_`j'
			if `nc' > 2 {
				qui gen tot_prob_lo = pp1lo_`j'
				qui gen tot_prob_hi = pp1hi_`j'

				forvalues c = 2/`nc' {
					if `c' == `nc' {
						qui gen pp`c'lo_`j' = 1-tot_prob_lo
						label var pp`c'lo_`j' "Predicted prob of `dv'=`c' at min `iv'"

						qui gen pp`c'hi_`j' = 1-tot_prob_hi
						label var pp`c'hi_`j' "Predicted prob of `dv'=`c' at max `iv' "							

						qui gen diff`c'_`j'=pp`c'hi_`j'-pp`c'lo_`j'
						label var diff`c'_`j' "Change in predicted prob of `dv'=`c' over range of `iv'"
					}
					else {
						qui gen est`c'lo_`j' = est1lo_`j'+(v_`th`c'_place'*`stub'`th`c'_place')
						qui gen pp`c'lo_`j'=1/(1+exp(-est`c'lo_`j'))-tot_prob_lo
						label var pp`c'lo_`j' "Predicted prob of `dv'=`c' at min `iv'"
						replace tot_prob_lo = tot_prob_lo + pp`c'lo_`j'
						
						qui gen est`c'hi_`j' = est1hi_`j'+(v_`th`c'_place'*`stub'`th`c'_place')
						qui gen pp`c'hi_`j'=1/(1+exp(-est`c'hi_`j'))-tot_prob_hi
						label var pp`c'hi_`j' "Predicted prob of `dv'=`c' at max `iv'"							
						replace tot_prob_hi = tot_prob_hi + pp`c'hi_`j'
						
						qui gen diff`c'_`j'=pp`c'hi_`j'-pp`c'lo_`j'
						label var diff`c'_`j' "Change in predicted prob of `dv'=`c' over range of `iv'"
					}
					sum pp`c'lo_`j' pp`c'hi_`j' diff`c'_`j'
				}
				drop tot_prob_lo tot_prob_hi
			}			
			drop est1lo_`j' est1hi_`j'
		}
	}
	log close `dv'_results
	view `dv'_results.log
		
	// Create data for graphs
	foreach a of varlist `iv_l1_interact' {
		foreach b of varlist `iv_l2_interact' {
			qui tab `a'
			if r(r) <= 5 {
				qui levelsof `a', local(`a'_vals)
			}
			else {
				qui sum `a'
				local min_`a' = round(r(min), .01)
				local max_`a' = round(r(max), .01)
				qui centile `a'
				local med_`a' = round(r(c_1), .01)
				local `a'_vals = "`min_`a'' `med_`a'' `max_`a''"
			}
			
			qui sum `b'
			local min_`b'=r(min)
			local max_`b'=r(max)
			
			qui levelsof `b', local(`b'_levels)
			local bl = wordcount("``b'_levels'")
			if `bl' <= 5 {						/*This will give bad results if categories are not evenly spaced*/
				local npts = `bl' - 1
			}
			else {
				local npts=15
			}
			local inc=(`max_`b''-`min_`b'')/`npts'

			local c = regexr("`iv_l2_interact'","`b'","")
			local c = trim("`c'")
			local total_others = wordcount("`c'")
			
			local f = wordcount("``a'_vals'")
			if `nc'==2 local ff = (`f'*3)+1		/*Will be position of varying values of l2 term in matrix*/
			else local ff = (`f'*`nc'*3)+1

			
			forvalues ii = 1/`f' {
				local fff`ii' = "e1`ii', l1`ii', u1`ii'"
				if `nc'>2 {
					forvalues c = 2/`nc' {
						local fff`ii' = "`fff`ii''"+", e`c'`ii', l`c'`ii', u`c'`ii'"
					}
				}						
			}	

			local f2 = "0"
			forvalues j = 1/`f' {
				if `nc'==2 local f2 = "`f2'"+",0,0,0"
				else {
					forvalues c = 1/`nc' {
						local f2 = "`f2'"+",0,0,0"
					}
				}
			}
			matrix foo2 = `f2'
			
			forvalues iter = 0/`npts' {	
				local `b'T = `min_`b''+(`inc'*`iter')		
				local ii = 1
				foreach i2 in ``a'_vals' {
					capture drop est_others
					gen sub_others = 0
					gen add_others = 0
					
					forvalues oth = 1/`total_others' {
						local other = word("`c'",`oth')
						qui replace sub_others = sub_others + v_``other'_place'*(`stub'``other'_place' ///
							+ v_``a'_place'*`stub'``a'X`other'_place')					
						qui replace add_others = add_others + v_``other'_place'*(`stub'``other'_place' ///
							+ `i2'*`stub'``a'X`other'_place')
					}					
					qui gen est`ii' = base  ///
						- (v_``b'_place'*`stub'``b'_place'+ v_``a'_place'*`stub'``a'_place' ///
							+v_``a'_place'*v_``b'_place'*`stub'``a'X`b'_place'+sub_others)  ///
						+ (``b'T'*`stub'``b'_place'+`i2'*`stub'``a'_place' ///
							+`i2'*``b'T'*`stub'``a'X`b'_place'+add_others)
					qui gen pp1`ii'=1/(1+exp(-est`ii'))
					qui sum pp1`ii'
					scalar e1`ii'=r(mean)
					scalar l1`ii'=e1`ii'-1.95*r(sd)
					scalar u1`ii'=e1`ii'+1.95*r(sd)
					if `nc'>2 {
						qui gen tot_prob = pp1`ii'						
						forvalues c = 2/`nc' {
							if `c' == `nc' qui gen pp`c'`ii' = 1-tot_prob
							else {
								qui gen pp`c'`ii'=1/(1+exp(-est`ii'-(v_`th`c'_place'*`stub'`th`c'_place')))-tot_prob
								qui replace tot_prob = tot_prob + pp`c'`ii'
							}
							qui sum pp`c'`ii'
							scalar e`c'`ii'=r(mean)
							scalar l`c'`ii'=e`c'`ii'-1.95*r(sd)
							scalar u`c'`ii'=e`c'`ii'+1.95*r(sd)
						}
						drop tot_prob
					}
					drop sub_others add_others 
					drop est`ii' pp*`ii'
					local ++ii
				}
				*save values to matrix	
				matrix foo = `fff1'
				forvalues j = 2/`f' {
					matrix foo = foo, `fff`j''
				}
				matrix foo = foo, ``b'T'
				matrix foo2 = foo2 \ foo
			}
			
			
			*convert matrix to data
			matrix pts`a'X`b'=foo2[2..(`npts'+2),1..`ff']
			svmat pts`a'X`b'
			
			// Graphs 
			local ff0 = `ff'-1
			egen pp_max = rowmax(pts`a'X`b'1-pts`a'X`b'`ff0')
			qui sum pp_max
			local pp_max = ceil(r(max)*4)/4
			if `pp_max' == .25 local step = .05
			else if `pp_max' == .5 local step = .1
			else local step = .25
						
			if `nc'==2 { /*Graphs for dichotomous outcomes*/
				// Graphs for each value of l1 variable (for later use in graph combine)
				forvalues g = 1/`f' {
					local ggg = `g'*3-2
					local gg = `g'*3-1
					local g1 = `g'*3
					local av = word("``a'_vals'",`g')

					if `bl' <= 5 {	/*For categorical l2 var, graph effect only at actual values*/
						twoway rspike pts`a'X`b'`gg' pts`a'X`b'`g1' pts`a'X`b'`ff', ///
							legend(off) xtitle("`lab_`b''") ylabel(0(.25)1) xlabel(0(1)`max_`b'') ///
							ytitle("Predicted Probabilities", size(medlarge) margin(r+5) ///
							width(25)) graphregion(lcolor(white) color(white)) ///
							|| scatter pts`a'X`b'`ggg' pts`a'X`b'`ff', ///
							msize(medsmall) saving("`b'X`a'`av'.gph", replace) nodraw
						local graphs = "`graphs'"+" `b'X`a'`av'.gph"		
					}
					else {
						twoway connect pts`a'X`b'`ggg'-pts`a'X`b'`g1' pts`a'X`b'`ff',  ///
							lp(solid dash dash) lc(black gs8 gs8) m(i i i) ///
							legend(off) xtitle("`b'") yscale(range(0 `pp_max')) ylabel(#4) ///
							ytitle("Predicted Probabilities", size(medlarge) margin(r+5) ///
							width(25)) graphregion(lcolor(white) color(white)) ///
							saving("`b'X`a'`av'.gph", replace) nodraw
						local graphs = "`graphs'"+" `b'X`a'`av'.gph"
					}
				}
				// Superimposed graph for dichotomous l1 vars
				if `f' == 2 {
					if `bl' == 2 { /*For dichotomous l2 var, graph effects only at actual values*/
						twoway rspike pts`a'X`b'5 pts`a'X`b'6 pts`a'X`b'7, ///
							lc(black)  legend(off) xtitle("`lab_`b''") ylabel(0(.25)1) xlabel(0(1)1) ///
							ytitle("Predicted Probabilities", size(medlarge) ///
							margin(r+5) width(25)) graphregion(lcolor(white) color(white)) ///
							|| rspike pts`a'X`b'2 pts`a'X`b'3 pts`a'X`b'7, lc(gs8) ///
							|| scatter pts`a'X`b'1 pts`a'X`b'7, msize(medsmall) mc(black) ///
							|| scatter pts`a'X`b'4 pts`a'X`b'7, msize(medsmall) mc(gs8)
						graph save "`b'X`a'.gph", replace
					}
					else if `bl' > 5 {
						twoway connect pts`a'X`b'4-pts`a'X`b'6 pts`a'X`b'1-pts`a'X`b'3 pts`a'X`b'7, ///
							lp(solid solid solid dash dash dash) ///
							lc(black gs8 gs8 black gs8 gs8) m(i i i i i i) ///
							legend(off) xtitle("`lab_`b''") ylabel(0(.25)1) ///
							ytitle("Predicted Probabilities", size(medlarge)  ///
							margin(r+5) width(25)) graphregion(lcolor(white) color(white))
						graph save "`b'X`a'.gph", replace
					}
					if "`c(os)'"=="MacOSX" graph export "`b'X`a'.pdf", replace
					else graph export "`b'X`a'.tif", replace
				}
			}
			else {	/*Graphs for ordinal outcomes*/
				// Graphs for each outcome category for each value of l1 variable (for later use in graph combine)
				forvalues g = 1/`f' {
					forvalues dv_cat = 1/`nc' {
						local lastg = `g'-1
						local g1 = `nc'*`lastg'*3+`dv_cat'*3
						local gg = `g1'-1
						local ggg = `g1'-2
						local av = word("``a'_vals'",`g')					
						
						local al: value label `a'
						local aa = word("``a'_vals'", `g')
						capture local al1: label `al' `aa'
						if `dv_cat' == 1 local ycolor = "black"
						else local ycolor = "white"
						
						local dvl: value label `dv'
						capture local dvl1: label `dvl' `dv_cat'
						if `g' == `f' local tcolor = "black"
						else local tcolor = "white"
						
						if `bl' <= 5 {	/*For categorical l2 var, graph effect only at actual values*/
							twoway rspike pts`a'X`b'`gg' pts`a'X`b'`g1' pts`a'X`b'`ff', ///
								legend(off) xtitle("") ylabel(0(.25)1) xlabel(`min_`b''(1)`max_`b'') ///
								ytitle("`al1'", size(medlarge) margin(r+5) ///
									width(25) color(`ycolor'))  ///
								graphregion(lcolor(white) color(white)) ///
								|| scatter pts`a'X`b'`ggg' pts`a'X`b'`ff', ///
								msize(medsmall) saving("`b'X`a'`av'c`dv_cat'.gph", replace) nodraw
							local graphs`g' = "`graphs`g''"+" `b'X`a'`av'c`dv_cat'.gph"
						}
						else {
							twoway connect pts`a'X`b'`ggg'-pts`a'X`b'`g1' pts`a'X`b'`ff',  ///
								lp(solid dash dash) lc(black gs8 gs8) m(i i i) ///
								legend(off) xtitle("") ///
								ylabel(0(`step')`pp_max') ///
								ytitle("`al1'", size(medlarge) margin(r+5) width(25)  ///
									color(`ycolor')) ///
								title("`dvl1'", size(medlarge) color(`tcolor')) ///
								graphregion(lcolor(white) color(white)) ///
								saving("`b'X`a'`av'c`dv_cat'.gph", replace) nodraw
							local graphs`g' = "`graphs`g''"+" `b'X`a'`av'c`dv_cat'.gph"
						}
					}
				}
			}
//			// Combine graphs (if necessary)
//			if `f' > 2 {
//// Need to make specific graph combine commands here for each column so they can be labeled
//// Maybe label top-row graphs in advance . . . size an issue?
////				forvalues g = 1/`f' {
////					graph combine `graphs`g'', rows(1)
//				graph combine `graphs5' `graphs4' `graphs3' `graphs2' `graphs1',  ///
//					cols(`nc') imargin(zero)
//				if "`c(os)'"=="MacOSX" graph export "`b'X`a'.pdf", replace
//				else graph export "`b'X`a'.tif", replace
//			}
		}
	}
	saveold "`dv'_res.dta", replace
	cd ..
}	


//Create data for results dotplot
clear
set obs `total_ivs'
gen ivno = _n
save coeffplot.dta, replace



foreach dv in `dvs' {
	use "`dv'/`dv'_res.dta", clear
	quietly gen i_var = ""
	quietly gen ivno = .
	quietly gen b_`dv' = .
	quietly gen se_`dv' = .

	foreach a of varlist `iv_l1_interact' {
		gen `stub'``a'_place'a = `stub'``a'_place'
		foreach b of varlist `iv_l2_interact' {
			qui centile `b'
			qui replace `stub'``a'_place'a = `stub'``a'_place'a + `stub'``a'X`b'_place'*r(c_1)
			qui sum `b'
			qui replace `a'X`b'=(`b'-r(mean))*`a'
		}
	}
	
	local i = 1
	foreach v in `ivs' {
		qui replace i_var = "`lab_`v''" in `i'
		qui replace ivno = `i' in `i'
		qui sum `v'
		if `r(max)' == 1 {
			qui levelsof `v'
			if "`r(levels)'"=="0 1" local sdX2 = 1
			else {
				qui sum `v'
				local sdX2 = 2*`r(sd)'
			}
		}
		else local sdX2 = 2*`r(sd)'
		
		if !regexm("`iv_l1_interact'","`v'") qui sum `stub'`i' 
		else qui sum `stub'`i'a
		
		qui replace b_`dv' = `r(mean)'*`sdX2' in `i'
		if `nc' > 2 & `i'>=2 & `i'<=`total_ivs' qui replace b_`dv' = b_`dv'*-1 in `i' /*HLM quirk*/
		qui replace se_`dv' =`r(sd)'*`sdX2' in `i'	
	
		local ++i
	}
	
	keep ivno i_var b_`dv' se_`dv'
	drop if b_`dv'==.
	sort ivno
	save "`dv'/`dv'_est.dta", replace
	use coeffplot.dta, clear
	sort ivno
	merge ivno using "`dv'/`dv'_est.dta", _merge(_`dv')
	save coeffplot.dta, replace
}

drop if i_var==""
gen order=_n

append using `dataset'.dta

local i = -100
foreach b of varlist `iv_l2_interact' {
	replace order=`i' if order==``b'_place'
	local ++i
	foreach a of varlist `iv_l1_interact' {
		replace order=`i' if order==``a'X`b'_place'
		local ++i
	}
}
if "`iv_l1_interact'"=="income" & regexm("`iv_l2_interact'","gini_net") {
	sum order if i_var == "Years of Education"
	replace order=`r(mean)'+.5 if i_var=="Income Quintile"
}
sort order
drop if i_var=="Cons"
keep `dvs' i_var b_* se_* order
keep if order~=.
save coeffplot0.dta, replace

//Create results dotplot
use coeffplot0.dta, clear
sencode i_var, gen(iv_order)
reshape long b_ se_, i(i_var) j(dv) string
foreach v in `dvs' {
	qui replace dv = "`lab_`v''" if dv=="`v'"
}
sencode dv, gen(dv_order)
sort dv_order iv_order
gen ll = b_-1.96*se_
gen ul = b_+1.96*se_

if `total_dvs' == 1 {
	local o_s = 0
	local s_b = 0
}
else if `total_dvs' == 2 {
	local o_s = -.15
	local s_b = .3
}
else {
	local o_s = -.2
	local s_b = .2
}

qui sum iv_order
local n_iv = r(max)
qui sum ul
local ul_max = r(max)
qui sum ll
local ll_min = r(min)

local last_ll = ll in 1
local last_ul = ul in 1

if (`ul_max'-`last_ul')>abs(`ll_min'-`last_ll') local p = 1
else local p = 11

if `total_dvs' == 1 local leg="legend(off)"
else local leg = "legend(cols(1) pos(`p') ring(0) size(small) bmargin(vsmall) region(margin(small)))"

local l1_b = `total_ivs' - .75
local l1_t = `total_ivs' - `n_l1' - .4
local l2_b = `total_ivs' - `n_l1' - .6
local l2_t = `total_ivs' - `n_l1' - `n_l2' - .4
local l3_b = `total_ivs' - `n_l1' - `n_l2' - .6
local l3_t = `total_ivs' - `n_l1' - `n_l2' - `n_l3' - .4

local x1 = `ul_max'
local x2 = 1.025*`ul_max'

eclplot b_ ll ul iv_order if dv_order<=3, horizontal ylabel(1(1)`n_iv')  ///
	supby(dv_order, offset(`o_s') spaceby(`s_b'))  ///
	`leg'  ///
	rplottype(rspike) estopts(msize(small)) estopts1(mcolor(black))  ///
	estopts2(mfcolor(white)) estopts3(mfcolor(gs8)) ciopts(lcolor(black) lw(thin)) ///
	ytitle("") xtitle("") ///
	graphregion(lcolor(white) color(white)) xline(0, lpattern(dot))
  
graph save coeffplot.gph, replace
if "`c(os)'"=="MacOSX" graph export "coeffplot.pdf", replace
else graph export "coeffplot.tif", replace

save coeffplot.dta, replace

