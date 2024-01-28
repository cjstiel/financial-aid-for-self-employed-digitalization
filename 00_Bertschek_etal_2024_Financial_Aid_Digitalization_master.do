/*******************************************************************************
		  Article: German financial state aid during the COVID-19 pandemic
				- higher impact among digitalized self-employed? 
							
published in: Entrepreneurship & Regional Development, 2024, 36(1-2), 76-97.
authors: Irene Bertschek, Joern Block, Alexander S. Kritikos, Caroline Stiel	
affiliations: ZEW Mannheim, Trier University, DIW Berlin	
				
********************************************************************************
													                 
																	 Do-File 00
	
							MASTER DO-FILE
						
		CONTENT:	Root file that manages the execution of all 
					subordinated do-files.
		
		OUTLINE:	PART 1:	Prepare work space
					PART 2: Run programs

					
--------------------------------------------------------------------------------
code author: Caroline Stiel (DIW Berlin)
version: 11-Feb-2021
--------------------------------------------------------------------------------

********************************************************************************
							PART 1: Work space preparation
*******************************************************************************/

clear all
set more off

*=========================================*
* 1.1 directories
*=========================================*

	global root /* insert your root directory*/

	global input			"$root\01_input"
	global scripts			"$root\02_scripts"
	global results			"$root\04_results"

*=========================================*
* 1.2 packages
*=========================================*


/*
//for frequencies
	ssc install fre
	
// package kmatch 
	ssc install kmatch 
	ssc install moremata
	ssc install kdens
	
//for coefplot
	net install gr0059_1.pkg, from(http://www.stata-journal.com/software/sj15-1/)
	
//for pstest
cap ssc install psmatch2
*/	

********************************************************************************
* 					PART 2: Run programs
********************************************************************************
/*

*===================================================*
*	2.1 Tranform data
*===================================================*

	do "$scripts/01_Bertschek_etal_2024_Financial_Aid_Digitalization_transformation.do"

*===================================================*
*	2.2 Generate variables 
*===================================================*

	do "$scripts/02_Bertschek_etal_2024_Financial_Aid_Digitalization_variables.do"	
	
*===================================================*
*	2.3 Descriptive statistics 
*===================================================*

	do "$scripts/03_Bertschek_etal_2024_Financial_Aid_Digitalization_descriptives.do"

*================================================================*
*	2.4 Propensity score matching and treatment effects analysis
*================================================================*

	do "$scripts/04_Bertschek_etal_2024_Financial_Aid_Digitalization_treatmenteffects.do"
*/

		
********************************************************************************
* 								END
********************************************************************************		
