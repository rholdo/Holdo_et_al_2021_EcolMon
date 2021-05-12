# Holdo_et_al_2021_EcolMon
Source code and data for analysis of Serengeti seedling demography

The Analysis.R file replicates all analyses and figures in the paper
The code draws from the following data files:

1- Seedling_final.csv: raw data for all seedling/gulliver measurements
Column definitions:
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

