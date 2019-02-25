import os
import sys
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import query 
import useful_funcs
import math
from math import log10, radians, pi,cos,sin
import plot

class lc:

    def __init__(self,save_dir):
        
        self.lc_dir = "lightcurves/"
        useful_funcs.create_dir(self.lc_dir)
        self.query    = query.query_DES()
        self.query_S82 = query.Stripe82()
        self.query_SDSS = query.query_SDSS()
        self.save_dir = self.lc_dir+save_dir
        useful_funcs.create_dir(self.save_dir)
        self.quasar_catalog_dir = "catalog/"
        self.quasar_catalog = "DR14+DR7+OzDES+Milliq_S1S2.txt"
        self.lc_info_file = "lc_info.csv"
        self.log = "log"
        self.band_list = ["g","r","i","z"]


    def load_quasar_catalog(self):

        # loading the quasar catalog
        # columns=[ra,dec,z,flag1,flag2,flag3,flag4,mag_psf_i,spread_model_i,
        #           spread_model_err_i]
        print("loading quasar catalog: "+self.quasar_catalog_dir+\
              self.quasar_catalog)
        with open(self.quasar_catalog_dir+self.quasar_catalog) as f:
            first_line = f.readline()
        num_columns = len(first_line.split(","))
        flag_number = num_columns-6
        dtype = [('ra', '<f8'), ('dec', '<f8'), ('z', '<f8')]+\
                [("flag_"+str(i),"|S15") for i in range(flag_number)]+\
                [('mag_psf_i', '<f8'), ('spread_model_i', '<f8'),\
                 ('spread_model_err_i', '<f8')]
        quasars = np.genfromtxt(self.quasar_catalog_dir+self.quasar_catalog,\
                  delimiter=",",dtype=dtype)

        quasars.sort(order="ra")
        return quasars

    def write_header(self,quasar_catalog):

        f = open(self.lc_dir+self.lc_info_file,"w")
        headers = ",".join(quasar_catalog.dtype.names)
        f.write("name,"+headers+",N_DES_g,N_DES_r,N_DES_i,N_DES_z")
        f.write(",N_SDSS_g,N_SDSS_r,N_SDSS_i,N_SDSS_z\n")
        f.close()

    def get_zeropoint(self,filename,tag,data_coadd,band):
        
        data = self.query.get_firstcut_objects_from_filename(filename,tag,star=True)
        mag_psf = -2.5*np.log10(data["flux_psf"])
        data_clean = data[(mag_psf< 99) & (mag_psf> -99)]
        if len(data) == 0: return None,None
        coor_single = SkyCoord(ra=data_clean["ra"]*u.degree,
                               dec=data_clean["dec"]*u.degree)
        coor_coadd = SkyCoord(ra=data_coadd["ra"]*u.degree,
                             dec=data_coadd["dec"]*u.degree)
        match_sources = coor_single.match_to_catalog_sky(coor_coadd)
        match_indices = match_sources[0][match_sources[1]<1.0*u.arcsec]
        matched_data1 = data_clean[match_sources[1]<1.0*u.arcsec]
        matched_data2 = data_coadd[match_indices]
        mag_psf = -2.5*np.log10(matched_data1["flux_psf"])
        mag_diff = matched_data2["mag_psf_"+band]-mag_psf
        mag_diff_clean = mag_diff[(mag_diff<40)&(mag_diff>20)]
        #plot.plot_magnitude_comparison(matched_data2["mag_psf_"+band][(mag_diff<40)&(mag_diff>20)],mag_diff_clean,"plot/","test")
        print "Number of reference stars: ",len(mag_diff_clean)
        if len(mag_diff_clean)>=3:
            zeropoint = np.median(mag_diff_clean)
            zeropoint_rms = np.std(mag_diff_clean)/np.sqrt(len(mag_diff_clean))
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
        clean_objects = matched_quasars[(matched_quasars["flux_psf"]>0) &\
                                   (matched_quasars["flux_psf"]<10**10) &\
                                   (matched_quasars["mjd_obs"]>20000)&\
                                   (matched_quasars["mjd_obs"]<90000)]
        return clean_objects


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
        for filename in filename_list_clean["filename"]:
            print filename
            data = self.query.get_firstcut_objects_from_filename(filename,tag,star=False)
            band = filename.split("_")[1]
            zeropoint,zeropoint_rms = self.get_zeropoint(filename,tag,coadd_data,band)
            matched_quasar = self.get_target_quasar(quasar["ra"],quasar["dec"],\
                                                    data)
            if (zeropoint is None) or (zeropoint_rms is None):
                useful_funcs.print_and_write(self.log,"No zeropoint !!")
            elif len(matched_quasar) == 0:
                useful_funcs.print_and_write(self.log,"Quasar not found !!")
            else:
                quasar_cal = self.calibrate_mag(matched_quasar,zeropoint,zeropoint_rms)
                final_list = np.append(final_list,quasar_cal)
        clean_objects = final_list[(final_list["flux_psf"]>0) &\
                                   (final_list["flux_psf"]<10**10) &\
                                   (final_list["mjd_obs"]>20000)&\
                                   (final_list["mjd_obs"]<90000)]
        return clean_objects

    def convert_flux_to_mag(self,total_quasars):

        total_quasars["flux_err_psf"] = total_quasars["flux_err_psf"]*1.09/\
                                        total_quasars["flux_psf"]
        total_quasars["flux_psf"] = 22.5-2.5*np.log10(total_quasars["flux_psf"])
        total_quasars["flux_err_auto"] = total_quasars["flux_err_auto"]*1.09/\
                                        total_quasars["flux_auto"]
        total_quasars["flux_auto"] = 22.5-\
                                    2.5*np.log10(total_quasars["flux_auto"])
        total_quasars.dtype.names = tuple([w.replace("flux","mag") for w \
                                          in total_quasars.dtype.names])

        return total_quasars


    def generate_SDSS_lightcurve(self,quasar):

        dtype_SDSS = [("mjd_obs",float),("mag_psf",float),\
                      ("mag_err_psf",float),("band","|S1")]
        matched_quasars = np.array([],dtype=dtype_SDSS)
        objects = self.query_S82.q(quasar["ra"],quasar["dec"])
        if objects is None: return matched_quasars
        for band in self.band_list:
            objects_band = np.array(map(tuple,[list(row) + [band] \
                           for row in objects[["mjd_"+band,"mag_psf_"+band,\
                           "mag_err_psf_"+band]]]),dtype=dtype_SDSS)
            clean_objects = objects_band[(objects_band["mag_psf"]<30) &\
                                         (objects_band["mag_psf"]>10) &\
                                         (objects_band["mjd_obs"]>10000)&\
                                         (objects_band["mjd_obs"]<80000)]
            matched_quasars = np.append(matched_quasars,clean_objects)

        return matched_quasars

class spectra:

    def __init__(self):

        self.query_SDSS = query.query_SDSS()
        self.dir_filter = "filters/"  # create filter info dir
        useful_funcs.create_dir(self.dir_filter)
        self.dir_spec = "spectra/" # create dir for spectra
        useful_funcs.create_dir(self.dir_spec)
        self.band_list = ["g","r","i","z"]
        self.get_DES_SDSS_bandpass()


    def get_SDSS_spectrum(self,ra,dec,dist=4):

        ra_lower = ra - cos(radians(dec))*0.5*dist/3600.0 
        ra_upper = ra + cos(radians(dec))*0.5*dist/3600.0
        dec_lower = dec - 0.5*dist/3600.0
        dec_upper = dec + 0.5*dist/3600.0
        column,data = self.query_SDSS.main(["-q","select survey, run2d, "+\
                      "plate,mjd, fiberID, specObjID from "+\
                              "SpecObj where class = 'QSO' and ra between "+\
                              " %f and %f and dec between %f and %f " % \
                              (ra_lower,ra_upper,dec_lower, dec_upper)])
        if data is None: return None
        spec_data = [row.split("," )for row in data.split("\n")[:-1]]
        survey,run2d,plateid,mjd,fiberid,specid = spec_data[0]

        if survey == "boss": survey="eboss"

        url_spectrum = "https://data.sdss.org/sas/dr14/%s/spectro/redux/%s/"+\
                       "spectra/%s/spec-%s-%s-%s.fits"
        variables = (survey,run2d,plateid.zfill(4),plateid.zfill(4),\
                    mjd,fiberid.zfill(4))
        name = useful_funcs.degtohexname(ra,dec)
        useful_funcs.download_file(url_spectrum % variables, \
                                   self.dir_spec+name+".fits")

    def get_DES_SDSS_bandpass(self):
        
        info = {"SDSS":{},"DES":{}}

        # dowload SDSS bandpass
        url_filter_SDSS = "http://classic.sdss.org/dr3/instruments/imager/"+\
                          "filters/%s.dat"
        for band in self.band_list:
            useful_funcs.download_file(url_filter_SDSS % (band),
                                       self.dir_filter+"SDSS_"+band+".dat")

        # dowload DES bandpass
        url_filter_DES ="http://www.ctio.noao.edu/noao/sites/default/files/DECam/"+\
                        "STD_BANDPASSES_DR1.dat"
        useful_funcs.download_file(url_filter_DES,self.dir_filter+"DES.dat")

        # load SDSS bandpass info
        for band in self.band_list:
            data_SDSS = np.genfromtxt(self.dir_filter+"SDSS_"+band+".dat")
            info["SDSS"].update({band:data_SDSS[:,[0,1]]})

        # load DES bandpass info
        band_num = {"g":1,"r":2,"i":3,"z":4}
        data_DES = np.genfromtxt(self.dir_filter+"DES.dat")
        for band in self.band_list:
            info["DES"].update({band:data_DES[:,[0,band_num[band]]]})
        self.filter_info = info

    def spectrum_to_mag(self,lamb_spec,flux_spec,lamb_trans,ratio_trans):

        c = 3*10**10*10**8 # A/s
        ratio_spec = np.interp(lamb_spec,lamb_trans,ratio_trans)
        STL = flux_spec*lamb_spec*ratio_spec
        T_L = ratio_spec/lamb_spec
        f_nu = 1./c*np.trapz(STL,x=lamb_spec)/np.trapz(T_L,x=lamb_spec)
        m_AB = -2.5*np.log10(f_nu)-48.6

        return m_AB

    def mag_SDSS_to_DES(self,name):

        mag_diff = {}
        if not os.path.isfile(self.dir_spec+name+".fits"):
            for band in self.band_list : mag_diff.update({band:0.0})
            return mag_diff

        hdul = fits.open(self.dir_spec+name+".fits")
        data =  hdul[1].data
        data_clean = data[(data["ivar"]>0) & (data["and_mask"] == 0)]

        lambda_q = 10**data_clean["loglam"]
        flux_q =  10**(-17)*data_clean["flux"]
        for band in self.band_list :
            trans_SDSS,trans_DES = self.filter_info["SDSS"][band],\
                                   self.filter_info["DES"][band]
            mag_SDSS =  self.spectrum_to_mag(lambda_q,flux_q,trans_SDSS[:,0],\
                                             trans_SDSS[:,1])
            mag_DES =  self.spectrum_to_mag(lambda_q,flux_q,trans_DES[:,0],\
                                            trans_DES[:,1])
            #print ("mag_SDSS: "+str(mag_SDSS)+" ; mag_DES: "+str(mag_DES))
            diff = mag_DES-mag_SDSS
            print (band+": m_DES - m_SDSS =" + str(diff))
            if math.isnan(diff): mag_diff.update({band:0.0})
            else: mag_diff.update({band:diff})

        return mag_diff

class lc_single:


    def __init__(self):

        self.query = query.query_DES()
        self.band_list = ["g","r","i","z"]


    def generate_DES_lightcurves(self,ra,dec,model):

        matched_quasars = np.array([],dtype=self.query.dtype_single_model)
        objects = self.query.get_nearby_single_epoch_objects_model(\
                  ra,dec,10,model)
        if objects is None: return matched_quasars
        coor_quasar = SkyCoord(ra=ra*u.degree,\
                               dec=dec*u.degree)
        coor_objects = SkyCoord(ra=objects["ra"]*u.degree,\
                                dec=objects["dec"]*u.degree)
        dist = coor_quasar.separation(coor_objects)
        matched_quasars = objects[dist<1.0*u.arcsec ]
        clean_objects = matched_quasars[(matched_quasars["flux"]>0) &\
                                   (matched_quasars["flux"]<10**10) &\
                                   (matched_quasars["mjd_obs"]>20000)&\
                                   (matched_quasars["mjd_obs"]<90000)]
        return clean_objects

    def convert_flux_to_mag_model(self,quasars):

        quasars["flux_err"] = quasars["flux_err"]*1.09/quasars["flux"]
        quasars["flux"] = 22.5-2.5*np.log10(quasars["flux"])
        quasars.dtype.names = tuple([w.replace("flux","mag") for w \
                                     in quasars.dtype.names])
        return quasars




