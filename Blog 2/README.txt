**********************
Reference information:
**********************


Dataset title: Data for the article "Symbiotic nitrogen fixation does not stimulate soil phosphatase
	activity under temperate and tropical trees," published in 2023 by Oecologia.

Dataset contributors: Duncan N. L. Menge, Emily A. Jager, Andrew W. Quebbeman, Amelia A. Wolf, Steven 
	S. Perakis, and Jennifer L. Funk

Version number 1.0.0
Date: February 8, 2023
Dataset DOI: 10.XXX

Suggested citations:
	
	Dataset citation: 
	
	Menge DNL, Jager EA, Quebbeman AW, Wolf AA, Perakis SS, Funk JL. 2023. Data for the article
	"Symbiotic nitrogen fixation does not stimulate soil phosphatase activity under temperate and 
	tropical trees", Dryad, Dataset, https://doi.org/10.XXX

	Corresponding publications:

	Jager EA, Quebbeman AW, Wolf AA, Perakis SS, Funk JL, Menge DNL. 2023. Symbiotic nitrogen 
	fixation does not stimulate soil phosphatase activity under temperate and tropical trees. 
	Oecologia, in press.

Contact information:
	Duncan Menge
	Department of Ecology, Evolution, and Environmental Biology, Columbia University
	dm2972@columbia.edu
	ORCID ID: 0000-0003-4736-9844

Funding sources:
	National Science Foundation DEB-1457650 (DNLM, SSP)
	National Science Foundation DEB-1457444 (JLF)

---------------------------------------------------------------------------------------------------------------------

***************************************
Dates and locations of data collection:
***************************************

	Dates: Field data were collected between May 2015 and October 2018
	Locations: Field data were collected from Waiakea, Hawaii; Volcano, Hawaii; Starker Forest, Oregon;
		and Black Rock Forest, New York. See publication for coordinates and more details.

---------------------------------------------------------------------------------------------------------------------

****************************
Data and code file overview:
****************************

Summary metrics:
	File count: 23 files
	File breakdown and formats:
		1 README (this file, .txt)
		12 data files (.csv)
		10 code files (.R)
	Total file size: 0.75 MB

File naming conventions: Data files have a .csv suffix. Code files have a .R suffix. Files that pertain only to 
	individual sites have "HI_W," "HI_V," "OR," or "NY" prefixes to refer to Waiakea, Volcano, Oregon, and New
	York, respectively. Additional details on individual files are below. 

File names:

1. Data files:
A.	NY_FX_Size_FoliarCNIsotope_Data.csv (66 rows, 166 columns)
	OR_FX_Size_FoliarCNIsotope_Data.csv (64 rows, 122 columns)
	HI_W_FX_Size_FoliarCNIsotope_Data.csv (96 rows, 111 columns)
	HI_V_FX_Size_FoliarCNIsotope_Data.csv (108 rows, 111 columns)

B.	NY_FX_Soil_Pase.csv (66 rows, 13 columns)
	OR_FX_Soil_Pase.csv (64 rows, 10 columns)
	HI_W_FX_Soil_Pase.csv (108 rows, 11 columns)
	HI_V_FX_Soil_Pase.csv (96 rows, 11 columns)

C.	NY_FX_SNF.csv (66 rows, 51 columns)
	OR_FX_SNF.csv (64 rows, 44 columns)
	HI_W_FX_SNF.csv (108 rows, 29 columns)
	HI_V_FX_SNF.csv (96 rows, 29 columns)

2. Code files:

A.	NY_growth_analysis.R
	OR_growth_analysis.R
	HI_W_growth_analysis.R
	HI_V_growth_analysis.R

B.	SoilPase_Fig_AGB_AGR.R
	SoilPase_Fig_Pase_x_TRT.R
	SoilPase_Fig_Pase_x_TRT_Ndfa.R
	SoilPase_Fig_Pase_x_TRT_SNF.R
	Fig_3-6_setup.R

C.	SoilPase_PowerAnalysis.R

---------------------------------------------------------------------------------------------------------------------

***************************
Data and code file details:
***************************

In this section we describe (1) the data files, (2) the code files, and (3) the relationship between 
the code and the data. 

*********************************
1. The data files are as follows:
*********************************

****
 A.
****

NY_FX_Size_FoliarCNIsotope_Data.csv (66 rows, 166 columns)
OR_FX_Size_FoliarCNIsotope_Data.csv (64 rows, 122 columns)
HI_W_FX_Size_FoliarCNIsotope_Data.csv (96 rows, 111 columns)
HI_V_FX_Size_FoliarCNIsotope_Data.csv (108 rows, 111 columns)

General description: 

These are the datafiles that contain the raw data for measurements of size and foliar 
chemistry (%N and N isotopes) through time in each of the sites (NY = New York, OR = Oregon, 
HI_W = Waiakea, HI_V = Volcano). In the phosphatase paper they were used to calculate growth 
for Fig. 1.

Note about variable names: Many columns are repeated across multiple dates with suffixes
"_01," "_02", etc. The suffixes correspond to the measurement date period within a site. 
In the variable descriptions below we will only list each type of column once, with "_XX"
to indicate that there are multiple columns with the same type of data.

Variables:
	Counter: Row counter. Allowable values: natural numbers.
	Location: Which state. Allowable values: HI, OR, NY.
	Experiment: Code used for a larger project. Allowable value: FX.
	Site: Site within state. In Hawaii this means the actual site (Waiakea or Volcano). 
		In Oregon it means the pair grouping. In New York it means one of two locations
		within Black Rock Forest.
		Allowable values: Waiakea, Volcano, Turtle (in NY), Ski (in NY), natural numbers (in OR).
	Block (Hawaii sites only): Block grouping within site. Allowable values: natural numbers.
	PID (Plant ID): Unique identifier for each individual tree. Allowable values: natural numbers.
	Paired PID: PID of the paired tree for this tree. In Oregon and New York there are 
		fixer - non-fixer pairs, so each plant has a paired PID. In Hawaii there were
		fixer - fixer - non-fixer triads, so only the fixer paired PIDs (the non-fixers)
		are listed.
		Allowable values: natural numbers, NA.
	Species: Four-letter species codes. Allowable values: ROPS, BENI, ALRU, PSME, GLSE, CAEQ,
		PSCA, ACKO, MOFA, DOVI.
		ROPS: Robinia pseudoacacia
		BENI: Betula nigra
		ALRU: Alnus rubra
		PSME: Pseudotsuga menziesii
		GLSE: Gliricidia sepium
		CAEQ: Casuarina equisetifolia
		PSCA: Psidium cattleianum
		ACKO: Acacia koa
		MOFA: Morella faya
		DOVI: Dodonaea viscosa
	Treatment: Fertilization treatment. Note the these are coded differently in the
		datafiles compared to the paper. Allowable values: LN, MN, HN, PHN.
		LN = C ("low N" = "control")
		MN = +10 ("medium N" = "+10 g N/m2/y")
		HN = +15 ("high N" = "+15 g N/m2/y")
		PHN = +15+P ("phosphorus and high N" = "+15 g N/m2/y + 15 g P/m2/y")
	Lat_Deg: Degrees latitude for the tree (NY) or site (Waiakea). Not listed for Volcano or OR.
		Allowable values: real numbers, NA.
	Lon_Deg: Degrees longitude for the tree (NY) or site (Waiakea). Not listed for Volcano or OR.
		Allowable values: real numbers, NA.
	Date_XX: Date of sampling for this sampling period (the "XX") and this tree.
		Allowable values: dates, in the format MM/DD/YY.
	Days_XX: Number of days since the _01 date. Will be 0 for the _01 dates and above 0 for the rest.
		Allowable values: whole numbers.
	Diam_mm_XX: Basal diameter (as close to the ground as possible) in mm for this sampling period.
		Allowable values: decimal to tenths of mm or hundredths of mm when averages of two perpendicular
		measurements were taken, NA.
	Height_cm_XX: Height above ground in cm for this sampling period.
		Allowable values: natural numbers except when converted from inches, then decimal to hundredths
		of cm, NA.
	DBH_XX: Diameter at 1.3 m above the ground (DBH) in mm for this sampling period.
		Allowable values: decimal to tenths or hundredths of mm, NA.
	DBH_stem2_XX: DBH of a secondary stem in mm in this sampling period (NY _08 only)
		Allowable values: decimal to tenths or hundredths of mm, NA.
	DBH_stem3_XX: DBH of a secondary stem in mm in this sampling period (NY _08 only)
		Allowable values: decimal to tenths or hundredths of mm, NA.
	Canopy_length_cm_XX: Horizontal distance of one canopy axis in cm in this sampling period.
		Allowable values: natural numbers, NA.
	Canopy_width_cm_XX: Horizontal distance of the axis 90 degrees from length in cm in this sampling period.
		Allowable values: natural numbers, NA.
	Ground_cover_percent_XX: Visual estimate of the % cover of understory plants in this tree's plot
		in this sampling period. Ranges from 0-100. Only taken in some sampling periods, and
		not used in the paper.
		Allowable values: whole numbers 0-100, NA.
	Notes_XX: Free-form column for observations taken in the field for this tree in this sampling period.
		Allowable values: text strings, blank.
	Dead_XX: Code for whether plants were clearly alive (0), clearly dead (1), or somewhere in between.
		Allowable values: real numbers between 0 and 1 inclusive.
	Use_growth_XX: Code for whether to include in growth analysis based on tree's presence and health. 
		Anything with Use_growth_XX greater than or equal to 0.5 was used for growth analyses.
		Allowable values: real numbers between 0 and 1 inclusive.
	foliar_N_mg_g_XX: Nitrogen content (mg/g) from field-collected leaves from this sampling period.
		Allowable values: positive real numbers, NA.
	foliar_C_mg_g_XX: Carbon content (mg/g) from field-collected leaves from this sampling period.
		Allowable values: positive real numbers, NA.
	foliar_15N_AP_XX: 15N content (atom % 15N) from field-collected leaves from this sampling period.
		Allowable values: positive real numbers, NA.
	foliar_d15N_XX: 15N content (delta15N) from field-collected leaves from this sampling period.
		Allowable values: real numbers, NA.
	foliar_d13C_XX: 13C content (delta13C) from field-collected leaves from this sampling period.
		Allowable values: real numbers, NA.
	seedpod_N_mg_g_XX: Nitrogen content (mg/g) from field-collected leaves from this sampling period.
		(Robinia in NY only)
		Allowable values: positive real numbers, NA.
	seedpod_C_mg_g_XX: Carbon content (mg/g) from field-collected leaves from this sampling period.
		(Robinia in NY only)
		Allowable values: positive real numbers, NA.
	seedpod_15N_AP_XX: 15N content (atom % 15N) from field-collected leaves from this sampling period.
		(Robinia in NY only)
		Allowable values: positive real numbers, NA.
	seedpod_d15N_XX: 15N content (delta15N) from field-collected leaves from this sampling period.
		(Robinia in NY only)
		Allowable values: real numbers, NA.
	seedpod_d13C_XX: 13C content (delta13C) from field-collected seedpods from this sampling period. 
		(Robinia in NY only)
		Allowable values: real numbers, NA.
	disk_d15N_XX: 15N content (delta15N) from cellulose disks (average of all disks in the plot)
		from this sampling period.
		Allowable values: real numbers, NA.
	Aboveground_harvest_date_XX: Date of aboveground harvest. There was only a single harvest at 
		each site, but the sampling period code differed between NY and HI, so we're still
		listing it as "XX" here. (NY and HI only)
		Allowable values: dates, in the format MM/DD/YY.
	Harvest_notes_XX: Free-form column for observations during the harvest. (NY and HI only)
		Allowable values: text strings, blank.
	Main_stem_g_dry_XX: Dry biomass of the harvested main stem material in g. (NY and HI only)
		Although thorns are listed separately, this number includes thorn mass.
		Allowable values: positive real numbers, NA.
	Twigs_g_dry_XX: Dry biomass of the harvested twig material in g. (NY and HI only)
		Although thorns are listed separately, this number includes thorn mass.
		Allowable values: positive real numbers, NA.
	Leaves_g_dry_XX: Dry biomass of the harvested leaf material in g. (NY and HI only)
		Although thorns are listed separately, this number includes thorn mass.
		Allowable values: positive real numbers, NA.
	Secondary_stem_g_dry_XX: Dry biomass of the harvested secondary stem material in g. (NY and HI only)
		Although thorns are listed separately, this number includes thorn mass.
		Allowable values: positive real numbers, NA.
	Thorns_mg_g_dry_mainstem_XX: Dry biomass of the harvested thorn material on main stems 
		in mg (not g). (Robinia in NY only)
		Allowable values: positive real numbers, NA.
	Thorns_mg_g_dry_secondarystem_XX: Dry biomass of the harvested thorn material on secondary stems 
		in mg (not g). (Robinia in NY only)
		Allowable values: positive real numbers, NA.
	Thorns_mg_g_dry_twigs_XX: Dry biomass of the harvested thorn material on twigs  
		in mg (not g). (Robinia in NY only)
		Allowable values: positive real numbers, NA.
	Fruit_g_dry_XX: Dry biomass of the harvested fruit material in g. (HI only)
		Allowable values: positive real numbers, NA.
	Dead_g_dry_XX: Dry biomass of the harvested dead aboveground material in g. (HI only)
		Allowable values: positive real numbers, NA.
	Aboveground_mass_g_dry_XX: Total aboveground harvested dry biomass in g. (NY and HI only)
		Allowable values: positive real numbers, NA.
	Aboveground_mass_g_wet_XX: Total aboveground harvested wet biomass in g. (HI only)
		Allowable values: positive real numbers, NA.
	Root_g_dry_XX: Total belowground harvested dry biomass in g. (NY and HI only, but only used NY in the paper)
		Allowable values: positive real numbers, NA.
	Nodule_drymass_mg_XX: Nodule dry biomass in mg in cores (not upscaled to the whole tree). (NY only)
		Allowable values: positive real numbers, NA.
	Num_soil_cores_for_nodule_sampling_XX: The number of cores used for sampling nodules. (NY only)
		Allowable values: natural numbers, NA.
	Coring_radius_for_nodule_sampling_m_XX: The radius in m to which cores were taken for 
		sampling nodules. (NY only)
		Allowable values: positive real numbers, NA.
	Use_allometry_XX: Code for whether plants were used for constructing allometric equations.
		Values 0.3 and higher were used. (NY and HI only)
		Allowable values: real numbers between 0 and 1 inclusive.
	Allometry_notes_XX: Free-form column for observations about allometry. (NY and HI only)
		Allowable values: text strings, blank.
	Biomass_g_dry_XX: Total dry biomass in g of harvested material. (NY only)
		Allowable values: positive real numbers, NA.
	AGB_est_kg_XX: Estimated aboveground biomass in kg for sampling period XX based on size 
		and allometric equations.
		Allowable values: positive real numbers, NA.
	BGB_est_kg_XX: Estimated belowground biomass in kg for sampling period XX based on size
		and allometric equations. (NY and OR only)
		Allowable values: positive real numbers, NA.
	Biomass_est_kg_XX: Estimated total biomass in kg for sampling period XX based on size
		and allometric equations.(NY and OR only)
		Allowable values: positive real numbers, NA.
	Mainstem_est_kg_XX: Estimated main stem biomass in kg for sampling period XX based on size
		and allometric equations.(OR only)
		Allowable values: positive real numbers, NA.
	Leaves_est_kg_XX: Estimated leaf biomass in kg for sampling period XX based on size
		and allometric equations.(OR only)
		Allowable values: positive real numbers, NA.

****
 B. 
****

NY_FX_Soil_Pase.csv (66 rows, 13 columns)
OR_FX_Soil_Pase.csv (64 rows, 10 columns)
HI_W_FX_Soil_Pase.csv (108 rows, 11 columns)
HI_V_FX_Soil_Pase.csv (96 rows, 11 columns)

General description: 

These are the datafiles that contain soil phosphatase activity data. There are two sampling period (2017 and 2018)
for NY, compared to one for OR and HI.

Variables:
	Counter: Row counter. Allowable values: natural numbers.
	Location: Which state. Allowable values: HI, OR, NY.
	Experiment: Code used for a larger project. Allowable value: FX.
	Site: Site within state. In Hawaii this means the actual site (Waiakea or Volcano). 
		In Oregon it means the pair grouping. In New York it means one of two locations
		within Black Rock Forest.
		Allowable values: Waiakea, Volcano, Turtle (in NY), Ski (in NY), natural numbers (in OR).
	Block (Hawaii sites only): Block grouping within site. Allowable values: natural numbers.
	PID (Plant ID): Unique identifier for each individual tree. Allowable values: natural numbers.
	Species: Four-letter species codes. Allowable values: ROPS, BENI, ALRU, PSME, GLSE, CAEQ,
		PSCA, ACKO, MOFA, DOVI.
		ROPS: Robinia pseudoacacia
		BENI: Betula nigra
		ALRU: Alnus rubra
		PSME: Pseudotsuga menziesii
		GLSE: Gliricidia sepium
		CAEQ: Casuarina equisetifolia
		PSCA: Psidium cattleianum
		ACKO: Acacia koa
		MOFA: Morella faya
		DOVI: Dodonaea viscosa
	Treatment: Fertilization treatment. Note the these are coded differently in the
		datafiles compared to the paper. Allowable values: LN, MN, HN, PHN.
		LN = C ("low N" = "control")
		MN = +10 ("medium N" = "+10 g N/m2/y")
		HN = +15 ("high N" = "+15 g N/m2/y")
		PHN = +15+P ("phosphorus and high N" = "+15 g N/m2/y + 15 g P/m2/y")
	Date_XX: Date of sampling for this sampling period (the "XX") and this tree.
		Allowable values: dates, in the format MM/DD/YY.
	Use_growth_XX: Code for whether to include in growth analysis based on tree's presence and health. 
		Anything with Use_growth_XX greater than or equal to 0.5 was used for growth analyses.
		Allowable values: real numbers between 0 and 1 inclusive.
	Soil_Pase_umol_g_hr_XX: Soil phosphatase activity in umol per g dry soil per hr.
		Allowable values: non-negative real numbers, NA.

****
 C.
****

NY_FX_SNF.csv (66 rows, 51 columns)
OR_FX_SNF.csv (64 rows, 44 columns)
HI_W_FX_SNF.csv (108 rows, 29 columns)
HI_V_FX_SNF.csv (96 rows, 29 columns)

General description: 

These files have the symbiotic nitrogen fixation values, along with some relevant biomass measurements 
or estimates and foliar N measurements.

Note about variable names: Many columns are repeated across multiple dates with suffixes
"_01," "_02", etc. The suffixes correspond to the measurement date period within a site.
The measurements dates for each period are listed in NY_FX_Size_FoliarCNIsotope_Data.csv,
HI_W_FX_Size_FoliarCNIsotope_Data.csv, and HI_V_FX_Size_FoliarCNIsotope_Data.csv.
In the variable descriptions below we will only list each type of column once, with "_XX"
to indicate that there are multiple columns with the same type of data.

Variables:
	Counter: Row counter. Allowable values: natural numbers.
	Location: Which state. Allowable values: HI, OR, NY.
	Block (Hawaii sites only): Block grouping within site. Allowable values: natural numbers.
	PID (Plant ID): Unique identifier for each individual tree. Allowable values: natural numbers.
	Species: Four-letter species codes. Allowable values: ROPS, BENI, ALRU, PSME, GLSE, CAEQ,
		PSCA, ACKO, MOFA, DOVI.
	Treatment: Fertilization treatment. Note the these are coded differently in the
		datafiles compared to the paper. Allowable values: LN, MN, HN, PHN.
		LN = C ("low N" = "control")
		MN = +10 ("medium N" = "+10 g N/m2/y")
		HN = +15 ("high N" = "+15 g N/m2/y")
		PHN = +15+P ("phosphorus and high N" = "+15 g N/m2/y + 15 g P/m2/y")
	Date_XX: Date of sampling for this sampling period (the "XX") and this tree.
		Allowable values: dates, in the format MM/DD/YY.
	Use_growth_XX: Code for whether to include in growth analysis based on tree's presence and health. 
		Anything with Use_growth_XX greater than or equal to 0.5 was used for growth analyses.
		Allowable values: real numbers between 0 and 1 inclusive.
	foliar_N_mg_g_XX: Nitrogen content (mg/g) from field-collected leaves from this sampling period.
		Allowable values: positive real numbers, NA.
	Aboveground_mass_g_dry_XX: Total aboveground harvested dry biomass in g. (NY and HI only)
		Allowable values: positive real numbers, NA.
	Biomass_g_dry_XX: Total dry biomass in g of harvested material. (NY only)
		Allowable values: positive real numbers, NA.
	AGB_est_kg_XX: Estimated aboveground biomass in kg for sampling period XX based on size 
		and allometric equations.
		Allowable values: positive real numbers, NA.
	Biomass_est_kg_XX: Estimated total biomass in kg for sampling period XX based on size
		and allometric equations. (NY only)
		Allowable values: positive real numbers, NA.
	Ndfa_u_XX: %Ndfa (percent of N derived from N fixation) for each tree for sampling period XX,
		calculated from foliar N isotopes. This is the "base case" value from Menge et al. (2023)
		Ecological Monographs.
		Allowable values: real numbers, NA. Negative numbers are allowable because of how the isotope
		calculation works. See Menge et al. (2023) for details.
	Biomass_g_N_XX: Estimated total N in total biomass at sampling period XX, calculated from N
		concentrations and biomass in tissue types. See Menge et al. (2023) for details. (NY and OR)
		Allowable values: non-negative real numbers, NA.
	AGB_g_N_XX: Estimated total N in aboveground biomass at sampling period XX, calculated from N
		concentrations and biomass in tissue types. See Menge et al. (2023) for details. (HI only)
		Allowable values: non-negative real numbers, NA.
	Nfix_g_N_u_XX: g N fixed per tree by sampling period XX. See Menge et al. (2023) for details.
		Allowable values: non-negative real numbers, NA.
	Nfix_g_N_yr_u_XX: g N fixed per tree per year by sampling period XX.
		See Menge et al. (2023) for details.
		Allowable values: non-negative real numbers, NA.
	Nfix_g_N_yr_u_XX_YY: g N fixed per year between sampling periods XX and YY.
		See Menge et al. (2023) for details.
		Allowable values: non-negative real numbers, NA.
	Nfix_g_N_kg_biomass_yr_u_XX: g N fixed per kg biomass per year by sampling period XX.
		See Menge et al. (2023) for details.
		Allowable values: non-negative real numbers, NA.
	Nfix_g_N_kg_biomass_yr_u_XX_YY: g N fixed per kg biomass per year between sampling periods XX and YY.
		See Menge et al. (2023) for details.
		Allowable values: non-negative real numbers, NA.
	Nfix_g_N_g_biomass_N_yr_u_XX: g N fixed per g biomass N per year by sampling period XX.
		See Menge et al. (2023) for details.
		Allowable values: non-negative real numbers, NA.
	Nfix_g_N_g_biomass_N_yr_u_XX_YY: g N fixed per g biomass N per year between sampling periods XX and YY.
		See Menge et al. (2023) for details.
		Allowable values: non-negative real numbers, NA.
	Nfix_g_N_m2_yr_u_XX: g N fixed per square meter of canopy area per year by sampling period XX.
		See Menge et al. (2023) for details.
		Allowable values: non-negative real numbers, NA.
	Nfix_g_N_m2_yr_u_XX_YY: g N fixed per square meter of canopy area per year between sampling periods XX and YY.
		See Menge et al. (2023) for details.
		Allowable values: non-negative real numbers, NA.
	Notes: Free-form column for comments about these data.
		Allowable values: text strings, blank.

*********************************
2. The code files are as follows:
*********************************

The code files themselves are commented, so this README will not go into detail about each one.
Rather, it will give an overview of all the files.

*************************
A. Growth analysis files.
*************************

**
i.
**

NY_growth_analysis.R
OR_growth_analysis.R
HI_W_growth_analysis.R
HI_V_growth_analysis.R

These files read in the _Size_...csv datafiles, calculate relative and absolute growth rates 
for each time period of interest, and run the stats to determine pairwise differences. The 
analysis code prints out the statistical results, and the comments in the code record the 
letters from the pairwise differences. I.e., # a b b b means that the control treatment is 
significantly different from the +10, +15, and +15+P treatments, all of which are not 
significantly different from each other.

Note that these calculate stats for many more metrics than are included in the paper: each 
annual increment and the cumulative growth from the initial year to each successive year in 
addition to the cumulative growth from the initial year to the final year.



****************************
 B. Figure files.
****************************

These files call the analysis files when needed and construct the figures or print out information
for tables. Their names are self-explanatory. See the table below in (3) for more details on how they
link to the data files. These files also run the basic statistical tests reported in the paper.

SoilPase_Fig_AGB_AGR.R
SoilPase_Fig_Pase_x_TRT.R
SoilPase_Fig_Pase_x_TRT_Ndfa.R
SoilPase_Fig_Pase_x_TRT_SNF.R
Fig_3-6_setup.R
	

*****************************
 C. Power analysis file
*****************************

This file runs the power analysis reported in the paper.

SoilPase_PowerAnalysis.R

******************************************************************************************
3. The data and analysis files match to the display items as shown in the following table.
******************************************************************************************

Display Description		Data file(s)				Analysis code file(s)

Table 1	Site characteristics	N/A					N/A
Table 2	ANOVA: Pase ~ sp, trt	XX_FX_Soil_Pase.csv			SoilPase_Fig_Pase_x_TRT.R
Table 3	ANOVA: Pase ~ Ndfa, trt	XX_FX_Soil_Pase.csv			SoilPase_Fig_Pase_x_TRT_Ndfa.R
				XX_FX_SNF.csv
Table 4	ANOVA: Pase ~ Nfix, trt	XX_FX_Soil_Pase.csv			SoilPase_Fig_Pase_x_TRT_SNF.R
				XX_FX_SNF.csv
Fig. 1	Aboveground biomass	XX_FX_Size_FoliarCNIsotope_Data.csv	XX_growth_analysis.R
	absolute growth rate						Fig_3-6_setup.R
									SoilPase_Fig_AGB_AGR.R
Fig. 2	Pase ~ sp, trt		XX_FX_Soil_Pase.csv			SoilPase_Fig_Pase_x_TRT.R
Fig. 3	Pase ~ Ndfa, trt	XX_FX_Soil_Pase.csv			SoilPase_Fig_Pase_x_TRT_Ndfa.R
				XX_FX_SNF.csv
Fig. 4	Pase ~ Nfix, trt	XX_FX_Soil_Pase.csv			SoilPase_Fig_Pase_x_TRT_SNF.R
				XX_FX_SNF.csv

Note: "XX" in the filename means there are separate data or code files for NY, OR, HI_W, and HI_V.


*************
END OF README
*************