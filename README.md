# Project - Supermassive Black Hole Binary Finder

## Introduction

The goal of this code is to search for peioridc quasar candidates using light curves from Dark Energy Survey (DES), Sloan Digital Sky Survey (SDSS), and other archival imaging surveys. 

The regions of interest are DES Supernova(SN) fields S1 and S2. We will complie the light curves for quasars in DES-SN S1 and S2, and then search for possible periodic signals.

## Before running any analysis

### set PATH and PYTHONPATH

Need to run them first every time, if you want to run the script in a new terminal. 
If you don't want to set the path every time, you can replace `pwd` with the absolute path and put the two lines into your bash file (e.g., ~/.bashrc.)

```
export PATH=$PATH:`pwd`/bin
export PYTHONPATH=$PYTHONPATH:`pwd`/python
```

### Install eups 

In order to connect to DES (proprietary data) database, we need to install eups([manual link](https://opensource.ncsa.illinois.edu/confluence/display/DESDM/EUPS+User%27s+Guide)), which is a DESDM's package management system.
You might need a password for svn.
```
install_eups.sh
```

Follow the instruction and you can choose to intall eups at any path you want.

### Setup some eups packages

Need to execute this command before running the scripts in a new terminal every time.
For despydb, you need ~/.desservice.ini in order to connect to DES database.
```
source setup.sourceme
```
You can import more packages, but don't forget to install them first("eups distrib install $package $version") and put "setup $packgename" into setup.sourceme .

## Generate the quasar catalog 

### Download raw quasar catalogs.

Go to your workspace or where you want to store your data and do the analysis.
```
download_raw_catalog.sh
```
This will download 1. million quasar catalog, 2. OzDES catalog(No longer used, we query DR1 directely), 3. SDSS DR7 quasar catalog and 4. SDSS DR14 quasar catalog. Then, the code will cross-match and selecte spectroscapically confirmed quasars in DES-SN S1 and S2 fields.

### Process catalog

```
generate_quasar_catalog
```
The catalogs are in the "catalog/".
"DR14+DR7+OzDES+Milliq.txt" includes all the spectroscopically confirmed quasars in 10 DES-SN fields.
"DR14+DR7+OzDES+Milliq_S1S2.txt" includes all spectroscopically confirmed quasars in S1 and S2.

## Generate quasar light curves

Once we have the quasar catalog, we then make the optical light curves from DES and SDSS.
```
generate_lightcurves
```
This command will query DES and SDSS strip 82 database and download the single-epoch information. 
We also correct the magnitude difference due to different telescope systems.(tiny but important, SDSS-> DES).

```
lightcurve_stat
```
This command will tell us how many quasars have enough imaging epochs and the properties like magnitudes and number of epochs.
It will also generate the file "lightcurves/lc_clean.csv", which will be used as the final parent quasar catalog for the light curve analysis.

## Analyzing the light curves
Before you run the analyzing program, make sure run the program in background. (Take looooong time!, ~3days x 45-core )
you should understand how many cores your server has.
```
nohup lightcurve_analysis >> log_pring &
```
You can check "log_period_analysis" to see how many jobs are finished.

## Find strong candidates
```
periodogram_stat
```
Check "statistics/strong_candidates.csv" to see the strong candidates.

## Autocorrelation function (ACF)

We use the z-transformed discrete correlation function (ZCDF) code to perform ACF ([The website](https://webhome.weizmann.ac.il/home/tal/zdcf2.html)).






