import os
import sys
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import query 
import useful_funcs
#import plot

class lc:

    def __init__(self,save_dir):
        
        self.lc_dir = "lightcurves/"
        useful_funcs.create_dir(self.lc_dir)
        self.query    = query.query_DES()
        self.query_SDSS = query.Stripe82()
        self.save_dir = self.lc_dir+save_dir
        useful_funcs.create_dir(self.save_dir)
        self.quasar_catalog_dir = "catalog/"
        self.quasar_catalog = "spec_quasars_S1S2.txt"
        self.lc_info_file = "lc_info.csv"
        self.log = "log"
        self.band_list = ["g","r","i","z"]
        self.dtype = [("ra",float),("dec", float),("flag","|S4"),\
                      ("z",float),("where","|S6")]


    def load_quasar_catalog(self):

        # loading the quasar catalog
        # columns=[ra,dec,spec_flag,z,where,mag_psf_r,spread_model_r]
        print("loading quasar catalog: "+self.quasar_catalog_dir+\
              self.quasar_catalog)
        quasars = np.genfromtxt(self.quasar_catalog_dir+self.quasar_catalog,\
                  delimiter=",",dtype=self.dtype+[("mag_psf_r",float),\
                  ("spread_model_r",float)])

        quasars.sort(order="ra")
        return quasars

    def write_header(self):

        f = open(self.lc_dir+output,"w")
        f.write("name,ra,dec,redshift,mag_r,flag,N_g,N_r,N_i,N_z\n")
        f.close()

    def get_zeropoint(self,data,data_coadd,band):

        mag_psf = -2.5*np.log10(data["flux_psf"])
        data_clean = data[(mag_psf< 80) & (mag_psf< 80)]
        if len(data) == 0: return None,None
        coor_single = SkyCoord(ra=data_clean["ra"]*u.degree,
                               dec=data_clean["dec"]*u.degree)
        coor_coadd = SkyCoord(ra=data_coadd["ra"]*u.degree,
                             dec=data_coadd["dec"]*u.degree)
        match_sources = coor_single.match_to_catalog_sky(coor_coadd)
        match_indices = match_sources[0][match_sources[1]<1.0*u.arcsec]
        matched_data1 = data_clean[match_sources[1]<1.0*u.arcsec]
        matched_data2 = data_coadd[match_indices]
        #plot.plot_magnitude_comparison(matched_data1["MAG_PSF"],matched_data2["MAG_PSF_"+band],self.save_dir,mjd_obs)
        mag_psf = -2.5*np.log10(matched_data1["flux_psf"])
        mag_diff = matched_data2["mag_psf_"+band]-mag_psf
        mag_diff_clean = mag_diff[(mag_diff<40)&(mag_diff>20)]
        print "Number of reference stars: ",len(mag_diff_clean)
        if len(mag_diff)>=3:
            zeropoint = np.median(mag_diff)
            zeropoint_rms = np.std(mag_diff)/np.sqrt(len(mag_diff))
            return zeropoint,zeropoint_rms
        else: return None,None

    def get_target_quasar(self,ra,dec,data):

        coor_quasar = SkyCoord(ra=ra*u.degree,
                               dec=dec*u.degree)
        coor_objects = SkyCoord(ra=data["ra"]*u.degree,
                                dec=data["dec"]*u.degree)
        dist = coor_quasar.separation(coor_objects)
        matched_quasar = data[dist<2.0*u.arcsec ]

        return matched_quasar

    def calibrate_mag(self,quasar,zeropoint,zeropoint_rms):

        quasar["flux_err_psf"] = np.sqrt((zeropoint_rms*quasar["flux_psf"])**2+\
                                 quasar["flux_err_psf"]**2)\
                                 *10**(-0.4*zeropoint+9)
        quasar["flux_psf"] = quasar["flux_psf"]*10**(-0.4*zeropoint+9)
        quasar["flux_err_auto"] = np.sqrt((zeropoint_rms*quasar["flux_auto"])\
                                  **2+quasar["flux_err_auto"]**2)\
                                  *10**(-0.4*zeropoint+9)
        quasar["flux_auto"] = quasar["flux_auto"]*10**(-0.4*zeropoint+9)
        return quasar


    def generate_finalcut_lightcurve(self,quasar):

        matched_quasars = np.array([],dtype=self.query.dtype_single)
        objects = self.query.get_nearby_single_epoch_objects(\
                  quasar["ra"],quasar["dec"],10)
        if objects is None: return matched_quasars
        coor_quasar = SkyCoord(ra=quasar["ra"]*u.degree,\
                               dec=quasar["dec"]*u.degree)
        coor_objects = SkyCoord(ra=objects["ra"]*u.degree,\
                                dec=objects["dec"]*u.degree)
        dist = coor_quasar.separation(coor_objects)
        matched_quasars = objects[dist<1.0*u.arcsec ]
        return matched_quasars

    def generate_firstcut_lightcurve(self,quasar,year):

        tag = year+"_FIRSTCUT"
        # get magnitude of reference sources from coadd
        final_list = np.array([],dtype=self.query.dtype_single)
        coadd_data = self.query.get_nearby_coadd_objects(quasar["ra"],\
                     quasar["dec"],10*60)
        if coadd_data is None: return final_list
        useful_funcs.print_and_write(self.log,"Nearby Coadd objects:"+\
                                     str(len(coadd_data)))
        # get list of images
        filename_list = self.query.get_filename_list_from_tag(quasar["ra"],\
                        quasar["dec"],tag)
        if filename_list is None: return final_list
        filename_list.sort(order="id")
        filename_list = filename_list[::-1]
        expnums,idx_clean = np.unique(filename_list["expnum"],return_index=True)
        filename_list_clean = filename_list[idx_clean]
        # get magnitude of quasar in each image
        N = 1
        final_list = np.array([],dtype=self.query.dtype_single)
        for filename in filename_list["filename"]:
            print filename
            data = self.query.get_firstcut_objects_from_filename(filename,tag)
            band = filename.split("_")[1]
            zeropoint,zeropoint_rms = self.get_zeropoint(data,coadd_data,band)
            matched_quasar = self.get_target_quasar(quasar["ra"],quasar["dec"],\
                                                    data)
            if (zeropoint is None) or (zeropoint_rms is None):
                useful_funcs.print_and_write(self.log,"No zeropoint !!")
            elif len(matched_quasar) == 0:
                useful_funcs.print_and_write(self.log,"Quasar not found !!")
            else:
                quasar_cal = self.calibrate_mag(matched_quasar,zeropoint,zeropoint_rms)
                final_list = np.append(final_list,quasar_cal)
        return final_list


    def generate_SDSS_lightcurve(self,quasar):

        dtype_SDSS = [("mjd_obs",float),("mag_psf",float),\
                      ("mag_err_psf",float),("band","|S1")]
        matched_quasars = np.array([],dtype=dtype_SDSS)
        objects = self.query_SDSS.q(quasar["ra"],quasar["dec"])
        if objects is None: return matched_quasars
        for band in self.band_list:
            objects_band = np.array(map(tuple,[list(row) + [band] \
                           for row in objects[["mjd_"+band,"mag_psf_"+band,\
                           "mag_err_psf_"+band]]]),dtype=dtype_SDSS)
            clean_objects = objects_band[(objects_band["mag_psf"]<30) &\
                                         (objects_band["mag_psf"]>10) ]
            matched_quasars = np.append(matched_quasars,clean_objects)

        return matched_quasars


