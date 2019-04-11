#! /bin/bash -f 

# install eups (packages manager)
# you might need password for svn
# for despydb, you need desservice.ini

eups_dir="EUPS_DESDM"

if [ ! -e $eups_dir ]
then
    mkdir $eups_dir
fi

bin_path=$(dirname $(which install_eups.sh))
echo export PATH=$bin_path$'/:$PATH' >>setup.sourceme
echo export PYTHONPATH="${bin_path/bin/python}"$'/:$PYTHONPATH\n'>>setup.sourceme

cd $eups_dir

curl -O http://desbuild.cosmology.illinois.edu/desdm_eupsinstall.py
python desdm_eupsinstall.py

source ./eups/desdm_eups_setup.sh

#eups distrib install python 2.7.9+2
eups distrib install setuptools 18.0+1
eups distrib install numpy 1.10.4+4
eups distrib install scipy 0.14.0+10
eups distrib install despydb 2.0.4+0
eups distrib install despymisc 1.0.4+2
eups distrib install astropy 1.1.2+6
eups distrib install matplotlib 1.5.3+2
eups distrib install pandas  0.15.2+5
eups distrib install sextractor 2.24.4+1
eups distrib install psfex 3.21.0+5
eups distrib install swarp 2.40.1+0

cd ..

echo $'setup setuptools\n setup numpy\nsetup scipy\nsetup despydb'>>setup.sourceme
echo $'setup astropy\nsetup despymisc\nsetup matplotlib'>>setup.sourceme
echo $'setup pandas\nsetup sextractor\nsetup psfex\nsetup swarp\n'>>setup.sourceme

source setup.sourceme

## install some package not in EUPS

software_dir="software"
if [ ! -e $software_dir ]
then
    mkdir $software_dir
fi
cd $software_dir

pythonpackage_dir="python_package"
if [ ! -e $pythonpackage_dir ]
then
    mkdir $pythonpackage_dir
fi

# use python 2.7.9+2 in order to use pip without install
#setup python 2.7.9+2
#python -m ensurepip --default-pip
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python get-pip.py 

export PYTHONPATH=$PYTHONPATH:`pwd`/$pythonpackage_dir/lib/python2.7/site-packages/

echo $'export PYTHONPATH=$PYTHONPATH:'`pwd`/$pythonpackage_dir/$'lib/python2.7/site-packages/' >> ../setup.sourceme
echo $'export PATH=$PATH:'`pwd`/$pythonpackage_dir/$'bin/' >> ../setup.sourceme
echo $'export PATH=$PATH:'`pwd`/$'bin/' >> ../setup.sourceme

# install astroquery by pip
pip install --prefix=`pwd`/$pythonpackage_dir astroquery

# install astroML
pip install --prefix=`pwd`/$pythonpackage_dir astroML==0.3
pip install --prefix=`pwd`/$pythonpackage_dir astroML_addons==0.2.2

# install emcee
pip install --prefix=`pwd`/$pythonpackage_dir emcee

# install gatspy
pip install --prefix=`pwd`/$pythonpackage_dir gatspy

# intall scikit-learn
pip install --prefix=`pwd`/$pythonpackage_dir -U scikit-learn

# install  corner for ploting posterior
pip install --prefix=`pwd`/$pythonpackage_dir corner

# install carma by brandonckelly 
# https://github.com/brandonckelly/carma_pack/blob/master/Linux_install.txt
'''
mkdir boost
wget http://sourceforge.net/projects/boost/files/boost/1.54.0/boost_1_54_0.tar.bz2
tar --bzip2 -xf boost_1_54_0.tar.bz2
cd boost_1_54_0
./bootstrap.sh --prefix=`pwd`/../boost
./b2 install
cd ..
echo "export BOOST_DIR=`pwd`/boost"
export BOOST_DIR=`pwd`/boost
'''

eups distrib install boost 1.55.0+0
setup boost
echo $'setup boost\n' >> ../setup.sourceme

# 
mkdir armadillo
eups distrib install cmake 2.8.12.2+1
setup cmake
wget http://sourceforge.net/projects/arma/files/armadillo-9.200.7.tar.xz
tar xf armadillo-9.200.7.tar.xz
cd armadillo-9.200.7
cmake -DCMAKE_INSTALL_PREFIX=`pwd`/../armadillo -DBOOST_ROOT=$BOOST_DIR/  .
make install
cd ..
rm -rf armadillo-9.200.7*
mkdir `pwd`/armadillo/lib
ln -s `pwd`/armadillo/lib64/libarmadillo.so `pwd`/armadillo/lib/libarmadillo.so
export ARMADILLO_DIR=`pwd`/armadillo
echo "export ARMADILLO_DIR=$ARMADILLO_DIR" >> ../setup.sourceme 
export LD_LIBRARY_PATH=${ARMADILLO_DIR}/lib:${LD_LIBRARY_PATH}
echo $'export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ARMADILLO_DIR}/lib:${ARMADILLO_DIR}/lib64' >> ../setup.sourceme
mkl_file=`find /usr/local/intel/ -name libmkl_rt.so | grep libmkl_rt.so | grep 64 | head -1`
MKL_DIR="${mkl_file%/*}"
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$MKL_DIR
echo "export LD_LIBRARY_PATH="'${LD_LIBRARY_PATH}'":$MKL_DIR" >> ../setup.sourceme


if [[ ! $(git) ]]; then 
    eups install git 2.9.5+0
    setup git
fi
git clone https://github.com/brandonckelly/carma_pack.git
cd carma_pack/src
python setup.py install
cd ../..

pip install --prefix=`pwd`/$pythonpackage_dir acor





# install javelin
mkdir javelin
wget https://bitbucket.org/nye17/javelin/downloads/javelin-0.33.tar.gz
tar -xvzf javelin-0.33.tar.gz
cd javelin-0.33
python setup.py config_fc --fcompiler=gnu95 install --prefix=`pwd`/../javelin
cd ..
rm -rf javelin-0.33*
echo $'export PYTHONPATH=$PYTHONPATH:'`pwd`/javelin/lib/python2.7/site-packages/ >> ../setup.sourceme


# install fpack
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

