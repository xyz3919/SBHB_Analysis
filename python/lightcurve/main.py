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
        self.quasar_catalog_dir = "catalog/"
        self.quasar_catalog = "spec_quasars_S1S2.txt"
        self.lc_info_file = "lc_info.csv"

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

    def create_dir(self,directory):

        if not os.path.exists(directory):
            os.makedirs(directory)

