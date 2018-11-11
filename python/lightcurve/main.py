import os
import sys
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from query import query
import plot

class lc:

    def __init__(self,save_dir):
        
        self.lc_dir = "lightcurves/"
        if not os.path.isdir(self.lc_dir): os.mkdir(self.lc_dir)
        self.query    = query()
        self.save_dir = self.lc_dir+save_dir
        if not os.path.isdir(self.save_dir): os.mkdir(self.save_dir)
        self.quasar_catalog_dir = "catalog/"
        self.quasar_catalog = "spec_quasars_S1S2.txt"
        self.lc_info_file = "lc_info.csv"
        self.log = "log"

        self.dtype = [("ra",float),("dec", float),("flag","|S4"),\
                      ("z",float),("where","|S6")]


    def load_quasar_catalog(self):

        # loading the quasar catalog
        # columns=[ra,dec,spec_flag,z,where,mag_psf_r,spread_model_r]
        quasars = np.genfromtxt(self.quasar_catalog_dir+self.quasar_catalog,\
                  delimiter=",",dtype=self.dtype+[("mag_psf_r",float),\
                  ("spread_model_r",float)])

        quasars = fits.open(name)[1].data
        quasars.sort(order="ra")
        return quasars

    def write_header(self):

        f = open(self.lc_dir+output,"w")
        f.write("name,ra,dec,redshift,mag_r,flag,N_g,N_r,N_i,N_z\n")
        f.close()

    def generate_lightcurve(self,quasar):

        band_list = ["g","r","i","z"]
        name = useful_funcs.degtohexname(quasar["ra"],quasar["dec"])
        f = open(save_catalog,"a")
        f.write(name+",")
        f.write(",".join(np.array(quasar).astype(str)))
        self.create_dir(self.save_dir+name)
        matched_quasars = self.get_single_epoch_data(ra,dec)
        for band in band_list:
            print "- Fitting ",name,band
            matched_quasars_in_band = matched_quasars[matched_quasars["BAND"]==band]
            if len(matched_quasars_in_band) == 0 :
                print "No data found in "+band+" !"
                f.write(",0")
            else:
                f.write(","+str(len(matched_quasars_in_band)))
                matched_quasars_in_band.to_csv(self.save_dir+name+"/"+band+".csv",\
                                        index=False,columns=["MJD_OBS","FLUX_PSF",\
                                        "FLUX_ERR_PSF","FLUX_AUTO","FLUX_ERR_AUTO"])
        f.write("\n")
        f.close()
        return name
