#! /bin/bash -f 

# install eups (packages manager)
# you might need password for svn
# for despydb, you need desservice.ini

eups_dir="EUPS_DESDM"

if [ ! -e $eups_dir ]
then
    mkdir $eups_dir
fi

cd $eups_dir

curl -O http://desbuild.cosmology.illinois.edu/desdm_eupsinstall.py
python desdm_eupsinstall.py

eups distrib install despydb 2.0.4+0
eups distrib install  despymisc 1.0.4+2
eups distrib install astropy 1.1.2+6
eups distrib install scipy 0.14.0+12
eups distrib install matplotlib 1.5.3+2

