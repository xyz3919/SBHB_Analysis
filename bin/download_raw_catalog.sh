#! /bin/bash -f

raw_catalog_dir="raw_catalog"
catalog_dir="catalog"

## creat catalog directory

if [ ! -e $raw_catalog_dir ] || [ ! -e $catalog_dir ]
then
    mkdir $raw_catalog_dir
    mkdir $catalog_dir
fi

## download quasar catalog

# million quasar catalog
if [ ! -e $raw_catalog_dir/milliquas.txt ]
then
    wget http://quasars.org/milliquas.zip
    unzip milliquas.zip
    mv milliquas.txt ${raw_catalog_dir}/
    rm milliquas.zip
fi

# OzDES catalog
if [ ! -e $raw_catalog_dir/ajaa5b8dt2_mrt.txt ]
then
    wget http://iopscience.iop.org/1538-3881/153/3/107/suppdata/ajaa5b8dt2_mrt.txt
    mv ajaa5b8dt2_mrt.txt ${raw_catalog_dir}/
fi

## extract ra,dec,z,(spectroscopically confirmed , from where) 

# Million quasar catalog
cut -c 1-11,13-23,78-82,52-55,84-89 --output-delimiter "," $raw_catalog_dir/milliquas.txt > $raw_catalog_dir/million_quasar.txt
# OzDES catalog
awk 'NR>87{printf("%s,%s,%s\n",$3,$4,$5)}'  $raw_catalog_dir/ajaa5b8dt2_mrt.txt  > $raw_catalog_dir/OzDES.txt

