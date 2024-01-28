/*******************************************************************************
		  Article: German financial state aid during the COVID-19 pandemic
				- higher impact among digitalized self-employed? 
							
published in: Entrepreneurship & Regional Development, 2024, 36(1-2), 76-97.
authors: Irene Bertschek, Joern Block, Alexander S. Kritikos, Caroline Stiel	
affiliations: ZEW Mannheim, Trier University, DIW Berlin	
				
********************************************************************************
													                 
																	 Do-File 01
	
						TRANSFORM VARIABLES
				
					CONTENT: Rename variables	


			
--------------------------------------------------------------------------------
code author: Maximilian Priem (DIWecon)
version: 06-Sept-2022
--------------------------------------------------------------------------------					 
	 
*******************************************************************************/
set more off

* define dates and time etc.
* --------------------------
local date=ltrim("$S_DATE")
local date=subinstr("`date'"," ","_",2)

* start logfile
* -------------
cap log close
log using "$results\log_Tranformation_`date'.log", replace

* load data set
* -------------
use "$input/Data_2020_DIW_UTrier_ZEW_VGSD_SelfEmployed_Germany_Covid19.dta", clear

********************************************************************************	
* 					PART 1: Rename variables
********************************************************************************

rename q0001 umfang_se
rename q0002 start_se

rename q0003_0001 rechtsform_einzelunternehmen
rename q0003_0002 rechtsform_genossenschaft
rename q0003_0003 rechtsform_gbr
rename q0003_0004 rechtsform_gmbh
rename q0003_0005 rechtsform_ug
rename q0003_0006 rechtsform_verein
rename q0003_0007 rechtsform_sonst_per
rename q0003_0008 rechtsform_sonst_kap

rename q0004 Branche

rename q0005 employees
rename q0006 employees_450

rename q0007_0001 Arbeitsort_1
rename q0007_0002 Arbeitsort_2
rename q0007_0003 Arbeitsort_3
rename q0007_0004 Arbeitsort_4
rename q0007_0005 Arbeitsort_5
rename q0007_0006 Arbeitsort_6
rename q0007_0007 Arbeitsort_7
rename q0007_0008 Arbeitsort_8
rename q0007_0009 Arbeitsort_9
rename q0007_0010 Arbeitsort_10
rename q0007_0011 Arbeitsort_11

rename q0008_0001 Einschraenkung_1
rename q0008_0002 Einschraenkung_2
rename q0008_0003 Einschraenkung_3
rename q0008_0004 Einschraenkung_4
rename q0008_0005 Einschraenkung_5
rename q0008_0006 Einschraenkung_6
rename q0008_0007 Einschraenkung_7
rename q0008_0008 Einschraenkung_8
rename q0008_0009 Einschraenkung_9
rename q0008_0010 Einschraenkung_10
rename q0008_0011 Einschraenkung_11
rename q0008_other Einschraenkung_12

rename q0009_0001 Storno_volle_Ent
rename q0009_0002 Storno_teilw_Ent
rename q0009_0003 Storno_ohne_Ent
rename q0009_0004 Zahlungsverzoeg
rename q0009_0005 Zahlungsausf
rename q0009_0006 Nichts

rename q0010 Umsatzeinbrueche
rename q0011 Dauer_finanz_Durststr
rename q0012 Hoehe_Umsatzrueckgang
rename q0013 mtl_Betriebsausgaben
rename q0014 mtl_finanz_Luecke
rename q0015 running_cost
rename q0016 financing_gap
rename q0017 Plan_AG2
rename q0017 Plan_other 
rename q0018 household_size
rename q0019 private_liq_assets
rename q0020 Dauer_Zahlungsfaehig_betriebl
rename q0021 solvency_private

rename q0022_0001 Arbeitsweise_1
rename q0022_0002 Arbeitsweise_2
rename q0022_0003 Arbeitsweise_3
rename q0022_0004 Arbeitsweise_4
rename q0022_0005 Arbeitsweise_5
rename q0022_0006 Arbeitsweise_6
rename q0022_0007 Arbeitsweise_7
rename q0022_0008 Arbeitsweise_8
rename q0022_other Arbeitsweise_9

rename q0023_0001 Aenderung_1
rename q0023_0002 Aenderung_2
rename q0023_0003 Aenderung_3
rename q0023_0004 Aenderung_4
rename q0024 Dauer_Aenderung

rename q0025_0001 Digitalisierungsgrad_1
rename q0025_0002 Digitalisierungsgrad_2
rename q0025_0003 Digitalisierungsgrad_3
rename q0026_0001 Digitalisierung_veraendert_1
rename q0026_0002 Digitalisierung_veraendert_2
rename q0026_0003 Digitalisierung_veraendert_3

rename q0027_0001 Schadensbegr_1
rename q0027_0002 Schadensbegr_2
rename q0027_0003 Schadensbegr_3
rename q0027_0004 Schadensbegr_4
rename q0027_0005 Schadensbegr_5
rename q0027_0006 Schadensbegr_6
rename q0027_0007 Schadensbegr_7
rename q0027_0008 Schadensbegr_8
rename q0027_0009 Schadensbegr_9
rename q0027_0010 Schadensbegr_10
rename q0027_0011 Schadensbegr_11
rename q0027_0012 Schadensbegr_12
rename q0027_0013 Schadensbegr_13
rename q0027_other Schadensbegr_14

rename q0028_0001 staatl_Unt_1
rename q0028_0002 staatl_Unt_2 
rename q0028_0003 staatl_Unt_3   
rename q0028_0004 staatl_Unt_4   
rename q0028_0005 staatl_Unt_5 
rename q0028_0006 staatl_Unt_6 

rename q0029_0001 hilfr_Unt_1   
rename q0029_0002 hilfr_Unt_2     
rename q0029_0003 hilfr_Unt_3     
rename q0029_0004 hilfr_Unt_4     
rename q0029_0005 hilfr_Unt_5    
rename q0029_0006 hilfr_Unt_6 
 
rename q0030 Soforthilfe
rename q0031 Keine_Soforthilfe
rename q0031_other Keine_Soforthilge_other
rename q0032_0001 Antragsdatum
rename q0033 Antrag_Status
rename q0034 Grund_Ablehnung
rename q0034_other Grund_Ablehnung_other
rename q0035 Hilfe_ausgezahlt
rename q0035_other Tage_Hilfe_ausgezahlt
rename q0036 Tage_Antragstellung_Bewilligung
rename q0037_0001 Hilfen_gellungen
rename q0038_0001 Zufriedenheit_Umsetzung_SHilfe
rename q0039 Berechnung_SHilge
rename q0040_0001 Vergleich_Angestellte
rename q0041 Verbesserung_SHilfe

rename q0042 Miete_Buero
rename q0043_0001 Vermieter_1
rename q0043_0002 Vermieter_2
rename q0043_0003 Vermieter_3
rename q0044 Mietausfall

rename q0045_0001 Jahresvorraussage_1
rename q0045_0002 Jahresvorraussage_2

rename q0046_0001 Erl_Buerokr_1
rename q0046_0002 Erl_Buerokr_2
rename q0046_0003 Erl_Buerokr_3
rename q0046_0004 Erl_Buerokr_4
rename q0046_0005 Erl_Buerokr_5
rename q0046_0006 Erl_Buerokr_6
rename q0046_0007 Erl_Buerokr_7
rename q0046_0008 Erl_Buerokr_8
rename q0046_0009 Erl_Buerokr_9
rename q0046_other Erl_Buerokr_10

rename q0047 age
rename q0048 gender
rename q0049 location
rename q0050 Education
rename q0051_0001 risk


*==============================================================================*
* 							SAVE and CLEAN
*==============================================================================*

*save data
save "$input/Dataset_renamed.dta", replace

cap log close
********************************************************************************
* 									End
********************************************************************************