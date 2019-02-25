# Project - Supermassive Black Hole Binary Finder

## Before running any analysis

### setup PATH and PYTHONPATH

Need to run them first everytime, when you want to run the script in a new terminal. 
If you don't want to setup the path evrytime. You can replace the `pwd` with absolute path and put those two lines into your ~/.bashrc.
```
export PATH=$PATH:`pwd`/bin
export PYTHONPATH=$PYTHONPATH:`pwd`/python
```

### Install eups 

In order to connect to DES database, we need to install eups([manual link](https://opensource.ncsa.illinois.edu/confluence/display/DESDM/EUPS+User%27s+Guide)), which is a DESDM's new package management system.
You might need a password for svn.
```
install_eups.sh
```

Follow the instruction and you can choose to intall eups at any path you want.

### Setup some eups packages

Need to execute this command before running the scripts in a new terminal.
For despydb, you need ~/.desservice.ini in order to connect to DES database.
```
source setup.sourceme
```
You can set up more packages, but don't forget to install them first("eups distrib install $package $version") and put "setup $packgename" into setup.sourceme .

## Generate the quasar catalog 

### Download raw quasar catalogs.

Go to your workspace or where you want to store you data and do the analysis.
```
download_raw_catalog.sh
```
This will download 1. million quasar catalog, 2. OzDES catalog(No longer used, we query DR1 directely), 3. SDSS DR7 quasar catalog and 4. SDSS DR14 quasar catalog.

### Process catalog

```
generate_quasar_catalog
```
The catalogs are in the "catalog/".
"DR14+DR7+OzDES+Milliq.txt" includes all the spectroscopically confirmed quasars in 10 DES-SN fields.
"DR14+DR7+OzDES+Milliq_S1S2.txt" includes all spectroscopically confirmed quasars in S1 and S2.

## Generate quasar light curves

Once we have the quasar catalog, we then need to make the optical light curves from DES and SDSS.
```
generate_lightcurves
```
This program will query DES and SDSS strip 82 database and download the single-epoch information. 
We also correct the magnitude difference due to different telescope systems.(tiny but important, SDSS-> DES).

```
lightcurve_stat
```
This program will tell us how many quasar has enough imaging epochs and some of thier properties like magnitude and number of epochs.
It will also generate the file "lightcurves/lc_clean.csv", which will be used as final quasar catalog for the light curve analysis.

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





