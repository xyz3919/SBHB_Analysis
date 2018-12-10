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
eups distrib install pandas  0.15.2+5
eups distrib install sextractor 2.24.4+1
eups distrib install psfex 3.21.0+5
eups distrib install swarp 2.40.1+0

# install astroquery by pip

python -m ensurepip --default-pip

pip_used=`which pip`

$pip_used install astroquery

# install astroML

$pip_used install astroML
$pip_used install astroML_addons
$pip_used install emcee

# install fpack

software_dir="software"

if [ ! -e $software_dir ]
then
    mkdir $software_dir
fi

cd $software_dir
mkdir bin

wget http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio3450.tar.gz

tar xzf cfitsio3450.tar.gz
rm  cfitsio3450.tar.gz

cd cfitsio

./configure --prefix=`pwd`
make
make install
make fpack
make funpack
cp fpack ../bin/
cp funpack ../bin/


