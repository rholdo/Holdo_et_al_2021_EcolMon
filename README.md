# Holdo_et_al_2021_EcolMon
Source code and data for analysis of Serengeti seedling demography

The Analysis.R file replicates all analyses and figures in the paper
The code draws from the following data files:

1- Seedling_final.csv: Raw data for all seedling/gulliver measurements. Columns are as follows:
Subplot: unique seedling ID (two-letter transect code followed by two-digit number 00-99)
Transect: Transect name
Species: Tree species code (first 3 letters of genus followed by first 3 letters of species)
X: Location in Eastings (Datum: WGS84, projection: Arc_1960_36S)
Y: Location in Northings (Datum: WGS84, projection: Arc_1960_36S)
Per: Survey period (month-year)
Date: Date in m/d/yyyy
Year: Year
Day: Day of study, from January 30, 2016 (day 0)
Basal_cm: Basal diameter in cm
Ht_m: Height in m
Damage: Damaged (1) or not (0). Any damage, including fire and herbivory.
Alive: Alive (1) or dead (0)
Found: Found (1) or not (0) during a given survey
Topkilled: Topkilled (1) or not (0)
Fire: Fire damage (0 or 1, except from May 2018 onwards, where fire damage is scored on a scale from 0 to 10)
Herbivory: Browsing damage (0 or 1, except from May 2018 onwards, where fire damage is scored on a scale from 0 to 10)
Grass: Grass biomass in g/m2

2- Light_cal.csv: Calibration data for light sensor. Columns are as follows:
PAR: Photsynthetically active radiation in umol m-2 s-1
R: Resistance in ohms

3- Transect_TC_data.csv: Tree cover data for each seedling and sensor location. Columns are as follows:
ID: ID for seedling (=Subplot in Seedling_final.csv file), sensor, or camera (not used in this study)
Site: Trnasect name (=Transect in Seedling_final.csv file)
Type: Sensor, seedling or camera
X: Location in Eastings (Datum: WGS84, projection: Arc_1960_36S)
Y: Location in Northings (Datum: WGS84, projection: Arc_1960_36S)
TC30: Tree cover proportion obtained from Google Earth image classification within a 30-m circle centered in X,Y

4- Transect_Tree_Density_Data.csv: belt transect adult (> 2 m height) tree data. Columns are as follows:
Date: Survey date
Site: Transect name (=Transect in Seedling_final.csv file)
Tree_ID: unique ID for a given tree
X: Location in Eastings (Datum: WGS84, projection: Arc_1960_36S)
Y: Location in Northings (Datum: WGS84, projection: Arc_1960_36S)
Species: Tree species code (first 3 letters of genus followed by first 3 letters of species)
DBH1-4: Diameter at breast height (cm) of all stems
Basal 1-3: Basal diameter of all stems
Can1: Diameter of major crown axis (m)
Can2: Diameter of minor crown axis (m)
Ht: Tree height (m)

5- Transect_Seedling_Density_Data.csv: belt transect adult (> 2 m height) tree data. Columns are as follows:
Date: Survey date
Site: Transect name (=Transect in Seedling_final.csv file)
ID: unique ID for a given seedling
X: Location in Eastings (Datum: WGS84, projection: Arc_1960_36S)
Y: Location in Northings (Datum: WGS84, projection: Arc_1960_36S)
Species: Tree species code (first 3 letters of genus followed by first 3 letters of species)
Ht: Seedling height (cm)
Diam: Seedling basal diameter (cm)
Can1: Diameter of major crown axis (cm)
Can2: Diameter of minor crown axis (cm)
Seedl_Respr: whether thought to be a true seedling (S) or a resprout (R)

6- Daily_cleaned_microclimate_data.csv: Daily microclimate data. Columns are as follows:
ID: Sensor ID
Day: Day of study, from January 30, 2016 (day 0)
VWCmean: Mean daily volumetric water content in m3 m-3
Tmean: Mean daily soil surface temperature in C
Tmax: Maximum daily soil surface temperature in C
Tmin: Minimum daily soil surface temperature in C
PARmean: mean daily Photosynthetically active radiation in umol m-2 s-1
PARmax: maximum daily Photosynthetically active radiation in umol m-2 s-1
