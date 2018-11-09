# DES-SN_quasar_catalog

# setup PATH and PYTHONPATH

Need to run them first everytime you run the script in a new terminal. If yu don't want to setup the path evrytime. You can replace the `pwd` with absolute path and put it into you ~/.bashrc.
```
export PATH=$PATH:`pwd`/bin
export PYTHONPATH=$PYTHONPATH:`pwd`/python
```


# Install eups 

In order to connect to DES database, we need to eups([manual like](https://opensource.ncsa.illinois.edu/confluence/display/DESDM/EUPS+User%27s+Guide)), which is a DESDM's new package management system.
You might need password for svn

```
install_eups.sh
```

Follow the instruction and you can choose to intall eups at any path you want.


# Setup some eups packages

Need to run them before .
For despydb, you need ~/.desservice.ini in order to run it.

```
source setup.sourceme
```

You can setup more packages, but don't forget to install them first("eups distrib install $package $version") and put "setup $packgename" into setup.sourceme .
# Generate the quasar catalog 

## Download the raw quasar catalog.

Go to your workingspace or where you want to store you data 

```
download_raw_catalog.sh
```

This will download 1. million quasar catalog and  2. OzDES catalog.

## Process catalog

```
generate_quasar_catalog
```

The catalogs are in the catalog/ 
"milliqua+OzDES_SN.txt" includes all the quasars can quasar candidates in 10 DES-SN fields
"spec_quasars_S1S2.txt" includes all spectroscopically confirmed quasars in S1 and S2.




