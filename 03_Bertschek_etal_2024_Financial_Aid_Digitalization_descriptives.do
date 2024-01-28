/*******************************************************************************
		  Article: German financial state aid during the COVID-19 pandemic
				- higher impact among digitalized self-employed? 
							
published in: Entrepreneurship & Regional Development, 2024, 36(1-2), 76-97.
authors: Irene Bertschek, Joern Block, Alexander S. Kritikos, Caroline Stiel	
affiliations: ZEW Mannheim, Trier University, DIW Berlin	
				
********************************************************************************
													                 
																	 Do-File 03
							
							DESCRIPTIVE STATISTICS
							
		CONTENT:	Descriptive statistics
		
		OUTLINE:	PART 1: Definitions
					PART 2: Sociodemographics
					PART 3: Digitalization
					PART 4: Emergency Aid Program
					
					
--------------------------------------------------------------------------------
code authors: Caroline Stiel (DIW Berlin)
version: 21-October-2022 (v09)
--------------------------------------------------------------------------------
	
*******************************************************************************/
		
/*******************************************************************************
* PART 1: Definitions
*******************************************************************************/

set more off

*====================================*
* 1.1 Logfile
*====================================*

* define dates and time etc.
* -------------------------
local date=ltrim("$S_DATE")
local date=subinstr("`date'"," ","_",2)

* start log file
* -------------- 
cap log close
log using "$results\log_Descriptives_`date'.log", replace



*================================*
* 1.2 Load data
*================================*

* load cleaned data set of unmatched sample (16,859 obs)
* ------------------------------------------------------
use "$input/Datensatz_Deskriptives.dta", clear


********************************************************************************
* PART 2: Sociodemographics
********************************************************************************

* ---------------------------------------------------------
* Section 4.1 and Table A1: Self-employed persons by region
* ---------------------------------------------------------


* location
* --------
tab location, mi

* age
* ----
tab age_cat, mi

* gender
* ----- 
ab gender, mi

* education
* ---------
tab edu_cat, mi


********************************************************************************
* PART 3: Digitalization
********************************************************************************


* -------------------------------------------------------------------------
* Table 1: Digitalization level among German self-employed persons in 2020
* ------------------------------------------------------------------------


* pre-crisis digitalization level
* ---------------------------------
tab cat_dig1_vor_corona, mi
tab cat_dig2_vor_corona, mi
tab cat_dig3_vor_corona, mi

* level of digitalization during crisis
* -------------------------------------
tab Digitalisierung_veraendert_1, mi
tab Digitalisierung_veraendert_2, mi
tab Digitalisierung_veraendert_3, mi



********************************************************************************
* PART 4: Emergency Aid Program
********************************************************************************

* ------------------
* Section 4.1 (Data)
* ------------------

* number of applicants
* --------------------
tab Soforthilfe, mi

* application status
* ------------------
tab Antrag_Status2, mi

* obs with payment received
* -------------------------
tab Hilfe_ausgezahlt, mi



********************************************************************************
* 							Clean and save
********************************************************************************
	
cap log close	

********************************************************************************
* 									End
********************************************************************************	