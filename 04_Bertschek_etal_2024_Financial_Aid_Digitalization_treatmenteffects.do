/*******************************************************************************
		  Article: German financial state aid during the COVID-19 pandemic
				- higher impact among digitalized self-employed? 
							
published in: Entrepreneurship & Regional Development, 2024, 36(1-2), 76-97.
authors: Irene Bertschek, Joern Block, Alexander S. Kritikos, Caroline Stiel	
affiliations: ZEW Mannheim, Trier University, DIW Berlin	
				
********************************************************************************
													                 
																	 Do-File 04
							
							TREATMENT EFFECTS ANALYSIS
							
		CONTENT:	Propensity Score Matching and Treatment Effects Analysis
		
		OUTLINE:	PART 1: Definitions
					PART 2: Main model
					PART 3: Heterogeneous Effects Analysis
					PART 4: Nearest Neighbor Matching (Robustness Checks)

					
--------------------------------------------------------------------------------
code author: Caroline Stiel (DIW Berlin)
version: 24-March-2022 (v11)
--------------------------------------------------------------------------------
	
*******************************************************************************/

* define dates and time etc.
* -------------------------- 
local date=ltrim("$S_DATE")
local date=subinstr("`date'"," ","_",2)

* start logfile
* -------------
cap log close
log using "$results\log_Matching_`date'.log", replace


/*******************************************************************************
						PART 1: Definitions
*******************************************************************************/

*================================*
* 1.1 covariates for PSM
*================================*

* main model
* ----------
global attitudes "ib2.risk_c"
global crisisdem "i.running_cost_2 ib1.Hoehe_Umsatzrueckgang_2 ib1.solvency_firm dig_vor_corona Durststrecke3" 
global sociodem "se_full_time i.gender ib3.age_cat ib2.edu_cat ib2.location_DE i.Plan_AG2b week_quest" 
global busidem "ib2.duration_se ib1.industry_nace solo_solo" 

* heterogeneity analysis
* ----------------------
global crisisdem_dig "i.running_cost_2 ib1.Hoehe_Umsatzrueckgang_2 ib1.solvency_firm Durststrecke3"


*======================================*
* 1.2 number of bootstrap replications
*======================================*

* bootstrap replications
* ----------------------
global reps 1999



*================================*
* 1.3 Set input and output files
*================================*

* load data
* ---------
use "$input\Datensatz_Variables_v04.dta", clear

* set output file
* ---------------
global outfile "EmergencyAid_Matching_Digitalization_1999reps.xlsx"



*=====================*
* 1.4 outcome variable
*=====================*

* outcome variable in main analysis
* ---------------------------------
* binary variable
global depvar survival_di2 


/*******************************************************************************
						PART 2: Main model
*******************************************************************************/

*==============================================================================*	
* 2.1 Main model without trimming
*==============================================================================*

* prepare output file
* -------------------
putexcel set "$results/$outfile", modify sheet("Kernel Matching")
	putexcel C2 = "main model"
	putexcel C3 = "ATE"
	putexcel D3 = "ATT"
	putexcel B4 = "coef."
	putexcel B5 = "SE"
	putexcel B6 = "p-value"	
	putexcel B7 = "bs replications"
	putexcel B8 = "N matched"
	putexcel B9 = "N out of common support"
	putexcel B10 = "N total"
	putexcel B11 = "min"
	putexcel B12 = "max"

	
* do PSM Kernel matching
* ----------------------
kmatch ps Masn_e2 $sociodem $attitudes $crisisdem $busidem   ///
($depvar) if sampleT==1, ate att  pscmd(probit) ///
vce(bootstrap, rep($reps)) bwidth(cv) kernel(epan) 

* save results as matrix
* ----------------------
mat define M1=r(table)

* extract ATE, ATT
* --------------
mat M2 = M1[1,1..2]
mat list M2

* extract standard errors
* -----------------------
mat M3 = M1[2,1..2]
mat list M3

* extract p-value
* ---------------
mat M4 = M1[4,1..2]
mat list M4

* save all results as .xls
* -------------------------
putexcel set "$results/$outfile", modify sheet("Kernel Matching")
putexcel C4 = mat(M2)
putexcel C5 = mat(M3)
putexcel C6 = mat(M4)

* number of bootstrap replications
* ---------------------------------
putexcel C7= `e(N_reps)'

* number of matched observations
* ------------------------------
mat list e(_N)
mat define M5 = e(_N)
mat M5 = M5[3,1]
putexcel C8 = mat(M5) 
             
* number of observations
* ----------------------
putexcel C10= `e(N)'	
	
	
	
*==============================================================================*	
* 2.2 Main model with trimming
*==============================================================================*	

* --------------------------------
* Table 2: ATT for the main sample
* --------------------------------

* clean
* -----
cap: drop _KM_* 
scalar drop _all
matrix drop _all

* prepare output file
* -------------------
putexcel set "$results/$outfile", modify sheet("Kernel Matching")
putexcel E2 = "main model with (min,max) trimming"
putexcel E3 = "ATE"
putexcel F3 = "ATT"
putexcel G2 = "main model with (min,.95) trimming"
putexcel G3 = "ATE"
putexcel H3 = "ATT"

*==============================================================================*	
* 2.2.1 Obtain trimming interval
*==============================================================================*

* run initial PSM to obtain boundaries for trimming (min/max)
* -----------------------------------------------------------
kmatch ps Masn_e2 $sociodem $attitudes $crisisdem $busidem  ///
($depvar) if sample==1 , ate att  pscmd(probit) ///
bwidth(cv) kernel(epan)  gen(_KM_)

* Propensity score distribution of treatment group
* ------------------------------------------------
sum _KM_ps if Masn_e2==1, detail 
scalar sc_min1 = r(min)
scalar sc_max1 = r(max)
		
* Propensity score distribution of control group
* ----------------------------------------------
sum _KM_ps if Masn_e2==0, detail 
scalar sc_min2 = r(min)
scalar sc_max2 = r(max)

* choose lower bound as max(min treat, min control)
* ------------------------------------------------
if ( sc_min1 < sc_min2){
	local GMps_min = sc_min2
	}
	else{
		local GMps_min = sc_min1
		}

* choose upper bound as min(max treat, max control)
* ------------------------------------------------
if ( sc_max1 < sc_max2){
	local GMps_max = sc_max1
	}
	else{
		local GMps_max = sc_max2
		}
		
* display and save trimming boundaries
* ------------------------------------
display `GMps_min'		
display `GMps_max'
		

*==============================================================================*	
* 2.2.2 Run model with trimmed PSM distribution
*==============================================================================*

* do trimmed kernel PSM matching with endogenous boundaries
* ---------------------------------------------------------		
kmatch ps Masn_e2 $sociodem $attitudes $crisisdem $busidem   ///
($depvar) if sample==1, ate att comsup(`GMps_min' `GMps_max') pscmd(probit) ///
vce(bootstrap, rep($reps)) bwidth(cv) kernel(epan)		

* loop through two different max boundaries: GMps_max and .95
* -----------------------------------------------------------
local b  E G  
foreach a of numlist  `GMps_max' 0.95  {
		gettoken left b: b 
		display "`left'" 
		display "`b'" 

		* do trimmed kernel PSM matching
		* -------------------------------	
		kmatch ps Masn_e2 $sociodem $attitudes $crisisdem $busidem  ///
		($depvar) if sample==1, ate att comsup(`GMps_min' `a') pscmd(probit) ///
		vce(bootstrap, rep($reps)) bwidth(cv) kernel(epan) 

		* save results as matrix
		* ----------------------
		mat define M1=r(table)

		* extract ATE, ATT
		* --------------
		mat M2 = M1[1,1..2]
		mat list M2

		* extract standard errors
		* -----------------------
		mat M3 = M1[2,1..2]
		mat list M3

		* extract p-value
		* ---------------
		mat M4 = M1[4,1..2]
		mat list M4
			
		* save all results as .xls
		* ------------------------ 
		putexcel set "$results/$outfile", modify sheet("Kernel Matching")
			
		putexcel `left'4= mat(M2)
		putexcel `left'5= mat(M3)
		putexcel `left'6= mat(M4)
		
		* number of bootstrap replications
		* ---------------------------------
		putexcel `left'7= `e(N_reps)'

		* number of matched observations
		* ------------------------------
		mat list e(_N)
		mat define M5 = e(_N)
		mat M5 = M5[3,1]
		putexcel `left'8 = mat(M5) 
		
		* number of observations trimmed (out of common support)
		* -------------------------------------------------------
		putexcel `left'9 = `e(N_outsup)'
        
		* number of observations
		* ----------------------
		putexcel `left'10= `e(N)'
		
		* lower and upper bound
		* ---------------------
		putexcel `left'11= `GMps_min'
		putexcel `left'12=`a'

		}

		
		
*==============================================================================*	
* 2.3 Matching quality (main model)
*==============================================================================*

*==============================================================================*
* 2.3.1 t-tests for trimmed model
*==============================================================================*

* -------------------------------------
* Section 5.1: Matching quality
* -------------------------------------

* prepare output file
* --------------------
putexcel set "$results/$outfile", modify sheet("Matching_Quality")
putexcel B2 = "Matching Quality: Kernel Matching (Min/Max Trimming)"
putexcel C3 = "Before matching"
putexcel D3 = "After matching"

* do again PSM kernel matching for trimmed model
* ----------------------------------------------
cap: drop _KM_*
kmatch ps Masn_e2 $sociodem $attitudes $crisisdem $busidem  ///
($depvar) if sample==1, ate att  pscmd(probit) ///
bwidth(cv) kernel(epan) comsup(`GMps_min' `GMps_max') gen(_KM_) //vce bootstrap nicht mit gen kompatibel 

* t-Test, standardised percentage bias
* ------------------------------------
pstest $sociodem $attitudes $crisisdem $busidem  ///
, mweight(_KM_mw) treated(Masn_e2) label both 
return list 

* save results
* ------------
putexcel set "$results/$outfile", modify sheet("Matching_Quality")
local B = r(Baft)
local B2 = r(Bbef)
local R = r(Raft)
local R2 = r(Rbef)
local M = r(medbiasaft)
local M2 = r(medbiasbef)
local P = r(r2bef)
local P2 = r(r2aft)

putexcel B17 = "Rubins B in %"
putexcel C17 = `B2'
putexcel D17 = `B'

putexcel B18 = "Rubins R"
putexcel C18 = `R2'
putexcel D18 = `R'

putexcel B14 = "Mean absolute standardized bias in %"
putexcel C14 = `M2'
putexcel D14 = `M'

putexcel B21 = "Pseudo-RSquared"
putexcel C21 = `P'
putexcel D21 = `P2'


*==============================================================================*
* 2.3.2 overlap
*==============================================================================*

* ------------------------
* Figure 3: Common support
* ------------------------


local ps_low `GMps_min' `GMps_min' 0
local count = 1

foreach a of numlist `GMps_max' 0.95 1 {
	cap: drop _KM_*
	gettoken ps ps_low:ps_low
	
	* do kernel PSM matching for trimmed model
	* ----------------------------------------
	kmatch ps Masn_e2 $sociodem $attitudes $crisisdem $busidem  ///
	($depvar) if sample==1, ate att comsup(`ps' `a') pscmd(probit) ///
	gen(_KM_) bwidth(cv) kernel(epan)				
	
	* draw histogram of ps density for treatment and control group
	* -------------------------------------------------------------
	cap: drop fx0
	cap: drop fx1
	cap: drop x
	cap: drop fx
	kdensity _KM_ps, nograph generate(x fx) 
	kdensity _KM_ps if Masn_e2==1, nograph generate(fx0) at(x)
	kdensity _KM_ps if Masn_e2==0, nograph generate(fx1) at(x)
	label var fx0 "treatment group"
	label var fx1 "control group"
	line fx0 fx1 x, color(black black) lpattern (2 dash) ///
	lwidth (thick thick) sort ytitle(Density, size(medlarge) margin(medium)) ///
	graphregion(color(white)) xtitle("Propensity score", size(medlarge) margin(medium)) ///
	legend(order(1 "treatment group" 2 "control group") size(medlarge)) ///
	name(common_support_`count', replace)
	graph export "$results\Overlap_`count'.png", replace
	
	display "------------------------------------------------------------------"
	display "count `count' shows the distribution for the upper trimming bound  `a'"
	display "------------------------------------------------------------------"
	display "the lower trimming bound is `ps'"
	display "------------------------------------------------------------------"
	local ++count
}



*==============================================================================*	
*						PART 3: Heterogeneous effects analysis
*==============================================================================*

*==============================================================================*	
* 3.1 Digitalization level pre-crisis
*==============================================================================*


*==============================================================================*	
* 3.1.1 average across areas
*==============================================================================*

* -----------------------------------------------------------------
* Table 2: ATT by digitalization level before and during the crisis
* ------------------------------------------------------------------

* 1: low
* 2: medium
* 3: high

* clean
* ------
cap: drop _KM_* 
scalar drop _all
matrix drop _all


* prepare output file	
* -------------------
putexcel set "$results/$outfile", modify sheet("digitalized_pre (average)")
putexcel C1 = "av. digitalization level before crisis (Min/Max Trimming)"
putexcel C2 = "low" 
putexcel E2 = "medium"
putexcel G2 = "high"
putexcel B4 = "coef."
putexcel C3 = "ATE"
putexcel D3 = "ATT"
putexcel E3 = "ATE"
putexcel F3 = "ATT"
putexcel G3 = "ATE"
putexcel H3 = "ATT"
putexcel B5 = "SE"
putexcel B6 = "p-value"	
putexcel B7 = "bs replications"
putexcel B8 = "N matched"
putexcel B9 = "N out of common support"
putexcel B10 = "N total"
putexcel B11 = "min"
putexcel B12 = "max"

		* ----------------------------------------------------------------------
		* First step: Kernel PSM to obtain trimming parameters for digital model
		* ----------------------------------------------------------------------  
		cap: drop _KM_*
		kmatch ps Masn_e2 $sociodem $attitudes $crisisdem_dig $busidem  ///
		(survival_di2) if sample==1 , ///
		ate att over(cat_dig_vor_corona)  pscmd(probit) bwidth(cv) kernel(epan)  gen(_KM_)

		* ps distribution treatment group
		* -------------------------------
		sum _KM_ps if Masn_e2==1, detail 
		scalar sc_min1 = r(min)
		scalar sc_max1 = r(max)
		
		* ps distribution control group
		* -----------------------------
		sum _KM_ps if Masn_e2==0, detail 
		scalar sc_min2 = r(min)
		scalar sc_max2 = r(max)
		
		* choose lower bound as max(min treat, min control)
		* ------------------------------------------------
		if ( sc_min1 < sc_min2){
			local Kps_min = sc_min2
		}
		else{
			local Kps_min = sc_min1
		}
		
		* choose upper bound as min(max treat, max control)
		* ------------------------------------------------
		if ( sc_max1 < sc_max2){
			local Kps_max = sc_max1
		}
		else{
			local Kps_max = sc_max2
		}
		
		* display trimming boundaries
		* ---------------------------
		display `Kps_min'		
		display `Kps_max'
		
		* --------------------------------------------------
		* Second step: kernel PSM for trimmed digital model
		* --------------------------------------------------
		kmatch ps Masn_e2 $sociodem $attitudes $crisisdem_dig $busidem  ///
		(survival_di2) if sample==1, ///
		ate att over(cat_dig_vor_corona) comsup(`Kps_min' `Kps_max') pscmd(probit) bwidth(cv) ///
		kernel(epan) vce(bootstrap,reps($reps))
		
			* save results as matrix
			* ----------------------
			mat define M1=r(table)
			
			* loop through all digitization categories
			* ----------------------------------------
			local cells C E G 
			forvalues i =1/3 {
			gettoken left cells:cells
			
			local i1 = `i'*2-1
			local i2 = `i'*2
			local i3 = `i'*3
			
			* extract ATE, ATT
			* -----------------
			mat M2 = M1[1,`i1'..`i2'] 
			mat list M2
			
			* extract standard errors
			* -----------------------
			mat M3 = M1[2,`i1'..`i2']
			mat list M3
			
			* extract p-value
			* ---------------
			mat M4 = M1[4,`i1'..`i2']
			mat list M4
			
			* save all results as .xls 
			* -----------------------
			putexcel set "$results/$outfile", modify sheet("digitalized_pre (average)")
			
			putexcel `left'4= mat(M2)
			putexcel `left'5= mat(M3)	
			putexcel `left'6= mat(M4)	
			
			* number of bootstrap replications
			* ---------------------------------
			putexcel `left'7= `e(N_reps)'

			* number of matched observations
			* ------------------------------
			mat list e(_N)
			mat define M5 = e(_N)
			mat M6 = M5[`i3',1]
			putexcel `left'8 = mat(M6)
			
			* number of unmatched observations
			* ---------------------------------
			mat M7 = M5[`i3',2]
			putexcel `left'9 = mat(M7)
             
			* number of observations trimmed (out of common support)
			* -------------------------------------------------------
			putexcel `left'10 = `e(N_outsup)'
        
			* number of observations
			* ----------------------
			*putexcel `left'10= `e(N)'
			mat M8 = M5[`i3',3]
			putexcel `left'11 = mat(M8)
		
		
			* lower and upper bound
			* ---------------------
			putexcel `left'12= `Kps_min'
			putexcel `left'13= `Kps_max'
			}

			* ----------------------------------------------
			* Test for significant differences between ATTs
			* ----------------------------------------------
			
			* low and medium
			* ---------------
			lincom [1]ATT - [2]ATT
			test [1]ATT = [2]ATT
			
			* medium and high
			* ---------------
			lincom [2]ATT - [3]ATT
			test [2]ATT = [3]ATT
			
			* low and high
			* ------------
			lincom [1]ATT - [3]ATT
			test [1]ATT = [3]ATT

			
*==============================================================================*	
* 3.1.2 by area
*==============================================================================*


* ------------------------------------------
* Tables 3 and 4: ATT by digitalization area
* ------------------------------------------


* digitalization level i
* -----------------------
* 1: low
* 2: medium
* 3: high


* area p
* ---------
* 1: products and services
* 2: internal processes
* 3: customer relations and distribution channels

* clean
* -----
cap: drop _KM_* 
scalar drop _all
matrix drop _all

forvalues p=1/3 {
if `p' == 1 local CATEGORY products
if `p' == 2 local CATEGORY process
if `p' == 3 local CATEGORY sales

* prepare output file	
* -------------------
putexcel set "$results/$outfile", modify sheet("digitalized_pre (`CATEGORY')")
putexcel C1 = "`CATEGORY' digitalization before crisis (Min/Max Trimming)"
putexcel C2 = "low" 
putexcel E2 = "medium"
putexcel G2 = "high"
putexcel B4 = "coef."
putexcel C3 = "ATE"
putexcel D3 = "ATT"
putexcel E3 = "ATE"
putexcel F3 = "ATT"
putexcel G3 = "ATE"
putexcel H3 = "ATT"
putexcel B5 = "SE"
putexcel B6 = "p-value"	
putexcel B7 = "bs replications"
putexcel B8 = "N matched"
putexcel B9 = "N out of common support"
putexcel B10 = "N total"
putexcel B11 = "min"
putexcel B12 = "max"

		* ----------------------------------------------------------------------
		* First step: Kernel PSM to obtain trimming parameters for digitization model
		* ----------------------------------------------------------------------  
		cap: drop _KM_*
		kmatch ps Masn_e2 $sociodem $attitudes $crisisdem_dig $busidem  ///
		(survival_di2) if sample==1, ///
		ate att over(cat_dig`p'_vor_corona) pscmd(probit) bwidth(cv) kernel(epan)  gen(_KM_)

		* ps distribution treatment group
		* -------------------------------
		sum _KM_ps if Masn_e2==1, detail 
		scalar sc_min1 = r(min)
		scalar sc_max1 = r(max)
		
		* ps distribution control group
		* ------------------------------
		sum _KM_ps if Masn_e2==0, detail 
		scalar sc_min2 = r(min)
		scalar sc_max2 = r(max)
		
		* choose lower bound as max(min treat, min control)
		* ------------------------------------------------
		if ( sc_min1 < sc_min2){
			local Kps_min = sc_min2
		}
		else{
			local Kps_min = sc_min1
		}
		
		* choose upper bound as min(max treat, max control)
		* ------------------------------------------------
		if ( sc_max1 < sc_max2){
			local Kps_max = sc_max1
		}
		else{
			local Kps_max = sc_max2
		}
		
		* display trimming boundaries
		* ---------------------------
		display `Kps_min'		
		display `Kps_max'
		
		* --------------------------------------------------
		* Second step: kernel PSM for trimmed digital model
		* --------------------------------------------------
		kmatch ps Masn_e2 $sociodem $attitudes $crisisdem_dig $busidem  ///
		(survival_di2) if sample==1, ///
		ate att over(cat_dig`p'_vor_corona) comsup(`Kps_min' `Kps_max') pscmd(probit) bwidth(cv) ///
		kernel(epan) vce(bootstrap,reps($reps))
		
			* save results as matrix
			* ----------------------
			mat define M1=r(table)
			
			* loop through all digitization categories
			* ----------------------------------------
			local cells C E G 
			forvalues i =1/3 {
			gettoken left cells:cells
			
			local i1 = `i'*2-1
			local i2 = `i'*2
			local i3 = `i'*3
			
			* extract ATE, ATT
			* -----------------
			mat M2 = M1[1,`i1'..`i2'] 
			mat list M2
			
			* extract standard errors
			* -----------------------
			mat M3 = M1[2,`i1'..`i2']
			mat list M3
			
			* extract p-value
			* ---------------
			mat M4 = M1[4,`i1'..`i2']
			mat list M4
			
			* save all results as .xls 
			* ------------------------
			putexcel set "$results/$outfile", modify sheet("digitalized_pre (`CATEGORY')")
			
			putexcel `left'4= mat(M2)
			putexcel `left'5= mat(M3)	
			putexcel `left'6= mat(M4)	
			
			* number of bootstrap replications
			* ---------------------------------
			putexcel `left'7= `e(N_reps)'

			* number of matched observations
			* ------------------------------
			mat list e(_N)
			mat define M5 = e(_N)
			mat M6 = M5[`i3',1]
			putexcel `left'8 = mat(M6)
			
			* number of unmatched observations
			* ---------------------------------
			mat M7 = M5[`i3',2]
			putexcel `left'9 = mat(M7)
             
			* number of observations trimmed (out of common support)
			* -------------------------------------------------------
			putexcel `left'10 = `e(N_outsup)'
        
			* number of observations
			* ----------------------
			*putexcel `left'10= `e(N)'
			mat M8 = M5[`i3',3]
			putexcel `left'11 = mat(M8)
		
		
			* lower and upper bound
			* ---------------------
			putexcel `left'12= `Kps_min'
			putexcel `left'13= `Kps_max'
			}

			* ----------------------------------------------
			* Test for significant differences between ATTs
			* ----------------------------------------------
			
			* low and medium
			* ---------------
			lincom [1]ATT - [2]ATT
			test [1]ATT = [2]ATT
			
			* medium and high
			* ---------------
			lincom [2]ATT - [3]ATT
			test [2]ATT = [3]ATT
			
			* low and high
			* ------------
			lincom [1]ATT - [3]ATT
			test [1]ATT = [3]ATT
}		

		
		
*==============================================================================*	
* 3.2 digitalization level during crisis
*==============================================================================*

*==============================================================================*	
* 3.2.1 average across areas
*==============================================================================*

* -----------------------------------------------------------------
* Table 2: ATT by digitalization level before and during the crisis
* ------------------------------------------------------------------

* changes
* -------
* 1: decreased
* 2: unchanged
* 3: increased

* clean
* -----
cap: drop _KM_* 
scalar drop _all
matrix drop _all


* prepare output file	
* -------------------
putexcel set "$results/$outfile", modify sheet("digital_post (average)")
putexcel C1 = "av. digitalization level (Min/Max Trimming)"
putexcel C2 = "decreased" 
putexcel E2 = "unchanged"
putexcel G2 = "increased"
putexcel B4 = "coef."
putexcel C3 = "ATE"
putexcel D3 = "ATT"
putexcel E3 = "ATE"
putexcel F3 = "ATT"
putexcel G3 = "ATE"
putexcel H3 = "ATT"
putexcel B5 = "SE"
putexcel B6 = "p-value"	
putexcel B7 = "bs replications"
putexcel B8 = "N matched"
putexcel B9 = "N out of common support"
putexcel B10 = "N total"
putexcel B11 = "min"
putexcel B12 = "max"

		* ----------------------------------------------------------------------
		* First step: Kernel PSM to obtain trimming parameters for digital model
		* ----------------------------------------------------------------------  
		cap: drop _KM_*
		kmatch ps Masn_e2 $sociodem $attitudes $crisisdem $busidem  ///
		(survival_di2) if sample==1, ///
		ate att  over(cat_dig_nach_corona) pscmd(probit) bwidth(cv) kernel(epan)  gen(_KM_)

		* ps distribution treatment group
		* -------------------------------
		sum _KM_ps if Masn_e2==1, detail 
		scalar sc_min1 = r(min)
		scalar sc_max1 = r(max)
		
		* ps distribution control group
		* -----------------------------
		sum _KM_ps if Masn_e2==0, detail 
		scalar sc_min2 = r(min)
		scalar sc_max2 = r(max)
		
		* choose lower bound as max(min treat, min control)
		* ------------------------------------------------
		if ( sc_min1 < sc_min2){
			local Kps_min = sc_min2
		}
		else{
			local Kps_min = sc_min1
		}
		
		* choose upper bound as min(max treat, max control)
		* ------------------------------------------------
		if ( sc_max1 < sc_max2){
			local Kps_max = sc_max1
		}
		else{
			local Kps_max = sc_max2
		}
		
		* display trimming boundaries
		* ---------------------------
		display `Kps_min'		
		display `Kps_max'
		
		* --------------------------------------------------
		* Second step: kernel PSM for trimmed digital model
		* --------------------------------------------------
		kmatch ps Masn_e2 $sociodem $attitudes $crisisdem $busidem  ///
		(survival_di2) if sample==1, ///
		ate att over(cat_dig_nach_corona) comsup(`Kps_min' `Kps_max') pscmd(probit) bwidth(cv) ///
		kernel(epan) vce(bootstrap,reps($reps))
		
			* save results as matrix
			* ----------------------
			mat define M1=r(table)
			
			* loop through all digitization categories
			* ----------------------------------------
			local cells C E G 
			forvalues i =1/3 {
			gettoken left cells:cells
			
			local i1 = `i'*2-1
			local i2 = `i'*2
			local i3 = `i'*3
			
			* extract ATE, ATT
			* -----------------
			mat M2 = M1[1,`i1'..`i2'] 
			mat list M2
			
			* extract standard errors
			* -----------------------
			mat M3 = M1[2,`i1'..`i2']
			mat list M3
			
			* extract p-value
			* ---------------
			mat M4 = M1[4,`i1'..`i2']
			mat list M4
			
			* save all results as .xls 
			* ------------------------
			putexcel set "$results/$outfile",  modify sheet("digital_post (average)")
			
			putexcel `left'4= mat(M2)
			putexcel `left'5= mat(M3)	
			putexcel `left'6= mat(M4)	
			
			* number of bootstrap replications
			* ---------------------------------
			putexcel `left'7= `e(N_reps)'

			* number of matched observations
			* ------------------------------
			mat list e(_N)
			mat define M5 = e(_N)
			mat M6 = M5[`i3',1]
			putexcel `left'8 = mat(M6)
			
			* number of unmatched observations
			* ---------------------------------
			mat M7 = M5[`i3',2]
			putexcel `left'9 = mat(M7)
             
			* number of observations trimmed (out of common support)
			* -------------------------------------------------------
			putexcel `left'10 = `e(N_outsup)'
        
			* number of observations
			* ----------------------
			*putexcel `left'10= `e(N)'
			mat M8 = M5[`i3',3]
			putexcel `left'11 = mat(M8)
		
		
			* lower and upper bound
			* ---------------------
			putexcel `left'12= `Kps_min'
			putexcel `left'13= `Kps_max'
			}

			* ----------------------------------------------
			* Test for significant differences between ATTs
			* ----------------------------------------------
			
			* low and medium
			* ---------------
			lincom [1]ATT - [2]ATT
			test [1]ATT = [2]ATT
			
			* medium and high
			* ---------------
			lincom [2]ATT - [3]ATT
			test [2]ATT = [3]ATT
			
			* low and high
			* ------------
			lincom [1]ATT - [3]ATT
			test [1]ATT = [3]ATT
		

*==============================================================================*	
* 3.2.2 by area
*==============================================================================*

* ------------------------------------------
* Tables 3 and 4: ATT by digitalization area
* ------------------------------------------


* digitalization level i
* -----------------------
* 1: low
* 2: medium
* 3: high


* area p
* ---------
* 1: products and services
* 2: internal processes
* 3: customer relations and distribution channels


* clean
* ------
cap: drop _KM_* 
scalar drop _all
matrix drop _all

forvalues p=1/3 {
if `p' == 1 local CATEGORY products
if `p' == 2 local CATEGORY process
if `p' == 3 local CATEGORY sales

* prepare output file	
* -------------------
putexcel set "$results/$outfile", modify sheet("digitalized_post (`CATEGORY')")
putexcel C1 = "`CATEGORY' digitalization during crisis (Min/Max Trimming)"
putexcel C2 = "decreased" 
putexcel E2 = "unchanged"
putexcel G2 = "increased"
putexcel B4 = "coef."
putexcel C3 = "ATE"
putexcel D3 = "ATT"
putexcel E3 = "ATE"
putexcel F3 = "ATT"
putexcel G3 = "ATE"
putexcel H3 = "ATT"
putexcel B5 = "SE"
putexcel B6 = "p-value"	
putexcel B7 = "bs replications"
putexcel B8 = "N matched"
putexcel B9 = "N out of common support"
putexcel B10 = "N total"
putexcel B11 = "min"
putexcel B12 = "max"

		* ----------------------------------------------------------------------
		* First step: Kernel PSM to obtain trimming parameters for digitization model
		* ----------------------------------------------------------------------  
		cap: drop _KM_*
		kmatch ps Masn_e2 $sociodem $attitudes $crisisdem $busidem  ///
		(survival_di2) if sample==1, ///
		ate att over(Digitalisierung_veraendert_`p') pscmd(probit) bwidth(cv) kernel(epan)  gen(_KM_)

		* ps distribution treatment group
		* -------------------------------
		sum _KM_ps if Masn_e2==1, detail 
		scalar sc_min1 = r(min)
		scalar sc_max1 = r(max)
		
		* ps distribution control group
		* -----------------------------
		sum _KM_ps if Masn_e2==0, detail 
		scalar sc_min2 = r(min)
		scalar sc_max2 = r(max)
		
		* choose lower bound as max(min treat, min control)
		* ------------------------------------------------
		if ( sc_min1 < sc_min2){
			local Kps_min = sc_min2
		}
		else{
			local Kps_min = sc_min1
		}
		
		* choose upper bound as min(max treat, max control)
		* ------------------------------------------------
		if ( sc_max1 < sc_max2){
			local Kps_max = sc_max1
		}
		else{
			local Kps_max = sc_max2
		}
		
		* display trimming boundaries
		* ---------------------------
		display `Kps_min'		
		display `Kps_max'
		
		* --------------------------------------------------
		* Second step: kernel PSM for trimmed digital model
		* --------------------------------------------------
		kmatch ps Masn_e2 $sociodem $attitudes $crisisdem $busidem  ///
		(survival_di2) if sample==1, ///
		ate att over(Digitalisierung_veraendert_`p') comsup(`Kps_min' `Kps_max') pscmd(probit) bwidth(cv) ///
		kernel(epan) vce(bootstrap,reps($reps))
		
			* save results as matrix
			* ----------------------
			mat define M1=r(table)
			
			* loop through all digitization categories
			* -------------------------------------
			local cells C E G 
			forvalues i =1/3 {
			gettoken left cells:cells
			
			local i1 = `i'*2-1
			local i2 = `i'*2
			local i3 = `i'*3
			
			* extract ATE, ATT
			* -----------------
			mat M2 = M1[1,`i1'..`i2'] 
			mat list M2
			
			* extract standard errors
			* -----------------------
			mat M3 = M1[2,`i1'..`i2']
			mat list M3
			
			* extract p-value
			* ---------------
			mat M4 = M1[4,`i1'..`i2']
			mat list M4
			
			* save all results as .xls 
			* ------------------------
			putexcel set "$results/$outfile",  modify sheet("digitalized_post (`CATEGORY')")
			
			putexcel `left'4= mat(M2)
			putexcel `left'5= mat(M3)	
			putexcel `left'6= mat(M4)	
			
			* number of bootstrap replications
			* ---------------------------------
			putexcel `left'7= `e(N_reps)'

			* number of matched observations
			* ------------------------------
			mat list e(_N)
			mat define M5 = e(_N)
			mat M6 = M5[`i3',1]
			putexcel `left'8 = mat(M6)
			
			* number of unmatched observations
			* ---------------------------------
			mat M7 = M5[`i3',2]
			putexcel `left'9 = mat(M7)
             
			* number of observations trimmed (out of common support)
			* -------------------------------------------------------
			putexcel `left'10 = `e(N_outsup)'
        
			* number of observations
			* ----------------------
			*putexcel `left'10= `e(N)'
			mat M8 = M5[`i3',3]
			putexcel `left'11 = mat(M8)
		
		
			* lower and upper bound
			* ---------------------
			putexcel `left'12= `Kps_min'
			putexcel `left'13= `Kps_max'
			}

			* ----------------------------------------------
			* Test for significant differences between ATTs
			* ----------------------------------------------
			
			* low and medium
			* ---------------
			lincom [1]ATT - [2]ATT
			test [1]ATT = [2]ATT
			
			* medium and high
			* ---------------
			lincom [2]ATT - [3]ATT
			test [2]ATT = [3]ATT
			
			* low and high
			* ------------
			lincom [1]ATT - [3]ATT
			test [1]ATT = [3]ATT		
}
	
	
	

*==============================================================================*	
* 3.3 regional analysis
*==============================================================================*

* -----------------------
* Table 5: ATT by region
* -----------------------


* clean
* -----
cap: drop _KM_* 
scalar drop _all
matrix drop _all


* prepare output file	
* --------------------
putexcel set "$results/$outfile", modify sheet("regions_Stadt_rest")
putexcel C1 = "regions: Stadtstaaten vs. rest (Min/Max Trimming)"
putexcel C2 = "other" 
putexcel E2 = "Stadtstaaten"
putexcel B4 = "coef."
putexcel C3 = "ATE"
putexcel D3 = "ATT"
putexcel E3 = "ATE"
putexcel F3 = "ATT"
putexcel B5 = "SE"
putexcel B6 = "p-value"	
putexcel B7 = "bs replications"
putexcel B8 = "N matched"
putexcel B9 = "N unmatched"
putexcel B10 = "N out of common support"
putexcel B11 = "N total"
putexcel B12 = "min"
putexcel B13 = "max"

* loop through all location categories
* -------------------------------------
		
		* ----------------------------------------------------------------------
		* First step: Kernel PSM to obtain trimming parameters for location model
		* ----------------------------------------------------------------------  
		cap: drop _KM_*
		kmatch ps Masn_e2 $sociodem $attitudes $crisisdem $busidem  ///
		(survival_di2) if sample==1, ///
		ate att over(location_dig_index_2) pscmd(probit) bwidth(cv) kernel(epan)  gen(_KM_)

		* ps distribution treatment group
		* -------------------------------
		sum _KM_ps if Masn_e2==1, detail 
		scalar sc_min1 = r(min)
		scalar sc_max1 = r(max)
		
		* ps distribution control group
		* -----------------------------
		sum _KM_ps if Masn_e2==0, detail 
		scalar sc_min2 = r(min)
		scalar sc_max2 = r(max)
		
		* choose lower bound as max(min treat, min control)
		* ------------------------------------------------
		if ( sc_min1 < sc_min2){
			local Kps_min = sc_min2
		}
		else{
			local Kps_min = sc_min1
		}
		
		* choose upper bound as min(max treat, max control)
		* ------------------------------------------------
		if ( sc_max1 < sc_max2){
			local Kps_max = sc_max1
		}
		else{
			local Kps_max = sc_max2
		}
		
		* display trimming boundaries
		* ---------------------------
		display `Kps_min'		
		display `Kps_max'
		
		* --------------------------------------------------
		* Second step: kernel PSM for trimmed digital model
		* --------------------------------------------------
		kmatch ps Masn_e2 $sociodem $attitudes $crisisdem $busidem  ///
		(survival_di2) if sample==1, ///
		ate att over(location_dig_index_2) comsup(`Kps_min' `Kps_max') pscmd(probit) bwidth(cv) ///
		kernel(epan) vce(bootstrap,reps($reps))
		
			* save results as matrix
			* ----------------------
			mat define M1=r(table)
							
			local cells C E
			
			forvalues i=1/2 {

			gettoken left cells:cells
			
			local i1 = `i'*2-1
			local i2 = `i'*2
			local i3 = `i'*3
			
			* extract ATE, ATT
			* -----------------
			mat M2 = M1[1,`i1'..`i2'] 
			mat list M2
			
			* extract standard errors
			* -----------------------
			mat M3 = M1[2,`i1'..`i2']
			mat list M3
			
			* extract p-value
			* ---------------
			mat M4 = M1[4,`i1'..`i2']
			mat list M4
			
			* save all results as .xls 
			* ------------------------
			putexcel set "$results/$outfile", modify sheet("regions_Stadt_rest")
			
			putexcel `left'4= mat(M2)
			putexcel `left'5= mat(M3)	
			putexcel `left'6= mat(M4)	
			
			* number of bootstrap replications
			* ---------------------------------
			putexcel `left'7= `e(N_reps)'

			* number of matched observations
			* ------------------------------
			mat list e(_N)
			mat define M5 = e(_N)
			mat M6 = M5[`i3',1]
			putexcel `left'8 = mat(M6)
			
			* number of unmatched observations
			* ---------------------------------
			mat M7 = M5[`i3',2]
			putexcel `left'9 = mat(M7)
             
			* number of observations trimmed (out of common support)
			* -------------------------------------------------------
			putexcel `left'10 = `e(N_outsup)'
        
			* number of observations
			* ----------------------
			*putexcel `left'10= `e(N)'
			mat M8 = M5[`i3',3]
			putexcel `left'11 = mat(M8)
		
		
			* lower and upper bound
			* ---------------------
			putexcel `left'12= `Kps_min'
			putexcel `left'13= `Kps_max'
			}

			* ----------------------------------------------
			* Test for significant differences between ATTs
			* ----------------------------------------------
			
			* other and Stadtstaaten
			* -----------------------
			lincom [1]ATT - [2]ATT
			test [1]ATT = [2]ATT

			
			
*==============================================================================*	
*				PART 4: Nearest Neighbour Matching (Robustness Checks)
*==============================================================================*


*==============================================================================	
* 4.1 pre-crisis digitalization level (average across areas)
*==============================================================================

* ----------------------------------------------------------------------
* Table 6: ATT with nearest neighbor-matching
* ----------------------------------------------------------------------


* levels
* ------
* 1: low
* 2: medium
* 3: high

* Prepare .xlsx spread sheet	
* --------------------------
putexcel set "$results/$outfile", modify sheet("Nearest Neighbour Matching PRE")
putexcel C1 = "Nearest Neighbour Matching (dig before crisis)"
putexcel C2 = "low"
putexcel D2 = "medium"
putexcel E2 = "high"
putexcel C3 = "ATT"
putexcel D3 = "ATT"
putexcel E3 = "ATT"
putexcel B4 = "coef."
putexcel B5 = "SE"
putexcel B6 = "p-value"	
putexcel B7 = "N total"

* Note: since bootstrapping does not provide robust standard errors for 
* NN-matching (Abadie and Imbens, 2008), we use the 'teffects' command with SE 
* calculated according to Abadie and Imbens (2012).



*==============================================================================*	
* 4.1.1 Low and high pre-crisis digitalization level
*==============================================================================*


 * estimate ATT in a loop for each digitalization level (low and high)
* --------------------------------------------------------------------
local cells C E 
foreach i in 1 3 {
	gettoken left cells:cells
    
scalar drop _all
matrix drop _all
cap noi drop OCS*

* obtain out-of sample information for number of NN-requirement
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig_vor_corona==`i'), ///
atet nneighbor(2) osample(OCS1)

* obtain out-of sample information for number of NN-requirement
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig_vor_corona==`i') & OCS1!=1, ///
atet nneighbor(2) osample(OCS2)

* run final estimation	 
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig_vor_corona==`i') & OCS1!=1& OCS2!=1, ///
atet nneighbor(2) 

			
		* save results as matrix
		* ----------------------
		mat define M1=r(table)
			
		* extract ATT
		* -----------------
		mat M2 = M1[1,1] 
		mat list M2
			
		* extract standard errors
		* -----------------------
		mat M3 = M1[2,1]
		mat list M3
			
		* extract p-value
		* ---------------
		mat M4 = M1[4,1]
		mat list M4
			
		* save all results as .xls 
		* ------------------------
		putexcel set "$results/$outfile", modify sheet("Nearest Neighbour Matching PRE")
		
		putexcel `left'4= mat(M2)
		putexcel `left'5= mat(M3)	
		putexcel `left'6= mat(M4)	

             
		* number of matched observations
		* ------------------------------
		putexcel `left'7= `e(N)'

}      

* save full dataset
* -----------------
save "$input/TEMP1.dta", replace


*==============================================================================*	
* 4.1.2 Medium pre-crisis digitalization level
*==============================================================================*


* use full dataset
* ----------------
use "$input/TEMP1.dta", clear

 * estimate ATT for digitalization level 'medium' separately 
* ----------------------------------------------------------
* reason: the sample of cat_dig_vor_corona==2 is too small and multicolinearity
* issues appear (perfect treatment predictions). Hence, the observations with
* perfect prediction probabilities must be dropped beforehand, since otherwise 
* the tseffects() command does not work and senseless output is produced.

* drop observations with perfect treatment probabilities
* -------------------------------------------------------
drop if location_DE==12 & sample==1 & cat_dig_vor_corona==2
drop if running_cost_2==9 & sample==1 & cat_dig_vor_corona==2

scalar drop _all
matrix drop _all
cap noi drop OCS*

* obtain out-of sample information for number of NN-requirement
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig_vor_corona==2), ///
atet nneighbor(2) osample(OCS1)

* obtain out-of sample information for number of NN-requirement
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig_vor_corona==2) & OCS1!=1, ///
atet nneighbor(2) osample(OCS2)

* run final estimation	 
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig_vor_corona==2) & OCS1!=1& OCS2!=1, ///
atet nneighbor(2) 

			
		* save results as matrix
		* ----------------------
		mat define M1=r(table)
			
		* extract ATT
		* -----------------
		mat M2 = M1[1,1] 
		mat list M2
			
		* extract standard errors
		* -----------------------
		mat M3 = M1[2,1]
		mat list M3
			
		* extract p-value
		* ---------------
		mat M4 = M1[4,1]
		mat list M4
			
		* save all results as .xls 
		* ------------------------
		putexcel set "$results/$outfile", modify sheet("Nearest Neighbour Matching PRE")
		
		putexcel D4= mat(M2)
		putexcel D5= mat(M3)	
		putexcel D6= mat(M4)	

             
		* number of matched observations
		* ------------------------------
		putexcel D7= `e(N)'

		
*==============================================================================	
* 4.2 pre-crisis digitalization level by area
*==============================================================================

* -------------------------------------------------------------
* Tables A2 and A3: ATT by value chain component (NN2 matching)
* --------------------------------------------------------------


*==============================================================================	
* 4.2.1 all areas
*==============================================================================

* use full dataset
* ----------------
use "$input/TEMP1.dta", clear

* digitalization level i
* -----------------------
* 1: low
* 2: medium
* 3: high


* area p
* ---------
* 1: products and services
* 2: internal processes
* 3: customer relations and distribution channels

* clean
* -----
cap: drop _KM_* 
scalar drop _all
matrix drop _all

forvalues p=1/3 {
if `p' == 1 local CATEGORY products
if `p' == 2 local CATEGORY process
if `p' == 3 local CATEGORY sales

* prepare output file	
* -------------------
putexcel set "$results/$outfile", modify sheet("NN2 digitalized_pre (`CATEGORY')")
putexcel C1 = "`CATEGORY' digitalization before crisis (NN2)"
putexcel C2 = "low" 
putexcel D2 = "medium"
putexcel E2 = "high"
putexcel B4 = "coef."
putexcel C3 = "ATT"
putexcel D3 = "ATT"
putexcel E3 = "ATT"
putexcel B5 = "SE"
putexcel B6 = "p-value"	
putexcel B7 = "N total"

* Note: since bootstrapping does not provide robust standard errors for 
* NN-matching (Abadie and Imbens, 2008), we use the 'teffects' command with SE 
* calculated according to Abadie and Imbens (2012).



* loop through all digitalization levels
* ------------------------------------
local cells C D E 
forvalues i =1/3 {
	gettoken left cells:cells
    
scalar drop _all
matrix drop _all
cap noi drop OCS*

* obtain out-of sample information for number of NN-requirement
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig`p'_vor_corona==`i'), ///
atet nneighbor(2) osample(OCS1)

* obtain out-of sample information for number of NN-requirement
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig`p'_vor_corona==`i') & OCS1!=1, ///
atet nneighbor(2) osample(OCS2)

* run final estimation	 
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig`p'_vor_corona==`i') & OCS1!=1& OCS2!=1, ///
atet nneighbor(2) 

			
		* save results as matrix
		* ----------------------
		mat define M1=r(table)
		
		* extract ATT			
		* -----------------
		mat M2 = M1[1,1] 
		mat list M2
			
		* extract standard errors
		* -----------------------
		mat M3 = M1[2,1]
		mat list M3
			
		* extract p-value
		* ---------------
		mat M4 = M1[4,1]
		mat list M4
			
		* save all results as .xls 
		* ------------------------
		putexcel set "$results/$outfile", modify sheet("NN2 digitalized_pre (`CATEGORY')")
		
		putexcel `left'4= mat(M2)
		putexcel `left'5= mat(M3)	
		putexcel `left'6= mat(M4)	

             
		* number of matched observations
		* ------------------------------
		putexcel `left'7= `e(N)'

}      
}

*==============================================================================	
* 4.2.2 Products & high, processes & medium, sales & medium
*==============================================================================

* Repeat estimation for those with perfect prediction probability errors

* reason: the sample of cat_dig_vor_corona==2 is too small and multicolinearity
* issues appear (perfect treatment predictions). Hence, the observations with
* perfect prediction probabilities must be dropped beforehand, since otherwise 
* the tseffects() command does not work and senseless output is produced.

* 1) products and high pre-crisis digitalization level
* ---------------------------------------------------
use "$input/TEMP1.dta", clear
cap noi drop OCS*

* drop observations with perfect treatment probabilities
* -------------------------------------------------------
drop if location_DE==14 & sample==1 & cat_dig1_vor_corona==3

* obtain out-of sample information for number of NN-requirement
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig1_vor_corona==3), ///
atet nneighbor(2) osample(OCS1)

* obtain out-of sample information for number of NN-requirement
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig1_vor_corona==3) & OCS1!=1, ///
atet nneighbor(2) osample(OCS2)

* run final estimation	 
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig1_vor_corona==3) & OCS1!=1& OCS2!=1, ///
atet nneighbor(2) 

			
		* Save results as matrix
		* ----------------------
		mat define M1=r(table)
			
		* extract ATT
		* ------------
		mat M2 = M1[1,1] 
		mat list M2
			
		* extract standard errors
		* -----------------------
		mat M3 = M1[2,1]
		mat list M3
			
		* extract p-value
		* ---------------
		mat M4 = M1[4,1]
		mat list M4
			
		* save all results as .xls 
		* ------------------------
		putexcel set "$results/$outfile", modify sheet("NN2 digitalized_pre (products)")
		
		putexcel E4= mat(M2)
		putexcel E5= mat(M3)	
		putexcel E6= mat(M4)	

             
		* number of matched observations
		* ------------------------------
		putexcel E7= `e(N)'
		
		
* 2) processes and medium pre-crisis digitalization level
* -------------------------------------------------------
use "$input/TEMP1.dta", clear
cap noi drop OCS*

* drop observations with perfect treatment probabilities
* -------------------------------------------------------
drop if industry_nace==3 & sample==1 & cat_dig2_vor_corona==2

* obtain out-of sample information for number of NN-requirement
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig2_vor_corona==2), ///
atet nneighbor(2) osample(OCS1)

* obtain out-of sample information for number of NN-requirement
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig2_vor_corona==2) & OCS1!=1, ///
atet nneighbor(2) osample(OCS2)

* run final estimation	 
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig2_vor_corona==2) & OCS1!=1& OCS2!=1, ///
atet nneighbor(2) 

			
		* save results as matrix
		* ----------------------
		mat define M1=r(table)
			
		* extract ATT
		* -----------------
		mat M2 = M1[1,1] 
		mat list M2
			
		* extract standard errors
		* -----------------------
		mat M3 = M1[2,1]
		mat list M3
			
		* extract p-value
		* ---------------
		mat M4 = M1[4,1]
		mat list M4
			
		* save all results as .xls 
		* ------------------------
		putexcel set "$results/$outfile", modify sheet("NN2 digitalized_pre (process)")
		
		putexcel D4= mat(M2)
		putexcel D5= mat(M3)	
		putexcel D6= mat(M4)	

             
		* number of matched observations
		* ------------------------------
		putexcel D7= `e(N)'
		

* 3) sales and medium pre-crisis digitalization level
* -------------------------------------------------------
use "$input/TEMP1.dta", clear
cap noi drop OCS*

* drop observations with perfect treatment probabilities
* -------------------------------------------------------
drop if location_DE==14 & sample==1 & cat_dig3_vor_corona==2

* obtain out-of sample information for number of NN-requirement
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig3_vor_corona==2), ///
atet nneighbor(2) osample(OCS1)

* obtain out-of sample information for number of NN-requirement
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig3_vor_corona==2) & OCS1!=1, ///
atet nneighbor(2) osample(OCS2)

* run final estimation	 
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem_dig $busidem , probit) if sample==1 & (cat_dig3_vor_corona==2) & OCS1!=1& OCS2!=1, ///
atet nneighbor(2) 

			
		* save results as matrix
		* ----------------------
		mat define M1=r(table)
			
		* extract ATT
		* -----------------
		mat M2 = M1[1,1] 
		mat list M2
			
		* extract standard errors
		* -----------------------
		mat M3 = M1[2,1]
		mat list M3
			
		* extract p-value
		* ---------------
		mat M4 = M1[4,1]
		mat list M4
			
		* save all results as .xls 
		* ------------------------
		putexcel set "$results/$outfile", modify sheet("NN2 digitalized_pre (sales)")
		
		putexcel D4= mat(M2)
		putexcel D5= mat(M3)	
		putexcel D6= mat(M4)	

             
		* number of matched observations
		* ------------------------------
		putexcel D7= `e(N)'





*==============================================================================*	
* 4.3  digitalization level during crisis (average across areas)
*==============================================================================

* ----------------------------------------------------------------------
* Table 6: ATT with nearest neighbor-matching
* ----------------------------------------------------------------------


* use full dataset
* ----------------
use "$input/TEMP1.dta", clear

* digitalization level
* --------------------
* 1: decreased
* 2: unchanged
* 3: increased

* Prepare .xlsx spread sheet	
* --------------------------
putexcel set "$results/$outfile", modify sheet("Nearest Neighbour Matching POST")
putexcel C1 = "Nearest Neighbour Matching (dig during crisis)"
putexcel C2 = "decreased"
putexcel D2 = "unchanged"
putexcel E2 = "increased"
putexcel C3 = "ATT"
putexcel D3 = "ATT"
putexcel E3 = "ATT"
putexcel B4 = "coef."
putexcel B5 = "SE"
putexcel B6 = "p-value"	
putexcel B7 = "N total"

* Note: since bootstrapping does not provide robust standard errors for 
* NN-matching (Abadie and Imbens, 2008), we use the 'teffects' command with SE 
* calculated according to Abadie and Imbens (2012).


*==============================================================================*	
* 4.3.1 Unchanged and increased digitalization levels
*==============================================================================*

* use full dataset
* ----------------
use "$input/TEMP1.dta", clear

 * estimate ATT in a loop for unchanged and increased digitalization level
* ------------------------------------------------------------------------

local cells D E 
forvalues i =2/3 {
	gettoken left cells:cells
    
scalar drop _all
matrix drop _all
cap noi drop OCS*

* obtain out-of sample information for caliper = 0.05
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem $busidem , probit) if sample==1 & (cat_dig_nach_corona==`i'), ///
atet nneighbor(2) osample(OCS1)

* obtain out-of sample information for number of NN-requirement
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem $busidem , probit) if sample==1 & (cat_dig_nach_corona==`i') & OCS1!=1, ///
atet nneighbor(2) osample(OCS2)

* run final estimation
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem $busidem , probit) if sample==1 & (cat_dig_nach_corona==`i') & OCS1!=1& OCS2!=1, ///
atet nneighbor(2)

			
		* save results as matrix
		* ----------------------
		mat define M1=r(table)
			
		* extract ATT
		* -----------------
		mat M2 = M1[1,1] 
		mat list M2
			
		* extract standard errors
		* -----------------------
		mat M3 = M1[2,1]
		mat list M3
			
		* extract p-value
		* ---------------
		mat M4 = M1[4,1]
		mat list M4
			
		* save all results as .xls 
		* ------------------------
		putexcel set "$results/$outfile", modify sheet("Nearest Neighbour Matching POST")
		
		putexcel `left'4= mat(M2)
		putexcel `left'5= mat(M3)	
		putexcel `left'6= mat(M4)	

 
		* number of matched observations
		* ------------------------------
		putexcel `left'7= `e(N)'

}      

*==============================================================================*	
* 4.3.2 Decreased digitalization level
*==============================================================================*

* use full dataset
* ----------------
use "$input/TEMP1.dta", clear


* estimate ATT for digitalization level 'decreased' separately 
* ----------------------------------------------------------
* reason: the sample of cat_dig_nach_corona==1 is too small and multicolinearity
* issues appear (perfect treatment predictions). Hence, the observations with
* perfect prediction probabilities must be dropped beforehand, since otherwise 
* the tseffects() command does not work and senseless output is produced.

* drop observations with perfect treatment probabilities
* -------------------------------------------------------
drop if (location_DE==4 | location_DE==12 | location_DE==16) & sample==1 & cat_dig_nach_corona==1
drop if (running_cost_2==9 | running_cost_2==11) & sample==1 & cat_dig_nach_corona==1

* join categories with too few observations (Hoehe Umsatzrückgang)
* ----------------------------------------------------------------
* join category: 'no decline' & 'decline up to 25%'
tab Hoehe_Umsatzrueckgang_2 if cat_dig_nach_corona==1  & sample==1, mi
recode Hoehe_Umsatzrueckgang_2 (2=1) (3=2) (4=3) (5=4) (6=5)
label define lbRevdecl 1 " max. Rückgang bis 25%" 2 "26%-50%" 3 "51%-75%" 4 "76%-99%" ///
5 "100% (Keine Einnahmen mehr)"
label value Hoehe_Umsatzrueckgang_2 lbRevdecl
tab Hoehe_Umsatzrueckgang_2 if cat_dig_nach_corona==1  & sample==1, mi

scalar drop _all
matrix drop _all
cap noi drop OCS*

* obtain out-of sample information for caliper = 0.05
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem $busidem , probit) if sample==1 & (cat_dig_nach_corona==1), ///
atet nneighbor(2) osample(OCS1)

* obtain out-of sample information for number of NN-requirement
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem $busidem , probit) if sample==1 & (cat_dig_nach_corona==1) & OCS1!=1, ///
atet nneighbor(2) osample(OCS2)

* run final estimation
cap noi teffects psmatch (survival_di2) (Masn_e2 $sociodem $attitudes ///
$crisisdem $busidem , probit) if sample==1 & (cat_dig_nach_corona==1) & OCS1!=1& OCS2!=1, ///
atet nneighbor(2)

			
		* save results as matrix
		* ----------------------
		mat define M1=r(table)
			
		* extract ATT
		* -----------------
		mat M2 = M1[1,1] 
		mat list M2
			
		* extract standard errors
		* -----------------------
		mat M3 = M1[2,1]
		mat list M3
			
		* extract p-value
		* ---------------
		mat M4 = M1[4,1]
		mat list M4
			
		* save all results as .xls 
		* ------------------------
		putexcel set "$results/$outfile", modify sheet("Nearest Neighbour Matching POST")
		
		putexcel C4= mat(M2)
		putexcel C5= mat(M3)	
		putexcel C6= mat(M4)	

 
		* number of matched observations
		* ------------------------------
		putexcel C7= `e(N)'



********************************************************************************
*    							Clean
********************************************************************************

* delete temp datasets	
* ---------------------	
cap noi erase "$input/TEMP1.dta"


********************************************************************************
*    							End
********************************************************************************
cap log close
