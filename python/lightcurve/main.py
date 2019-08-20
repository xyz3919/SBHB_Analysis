import os
import numpy as np
from scipy.signal import medfilt
from astropy.io import fits,ascii
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib
matplotlib.use('Agg')
import query 
import query_PS # get PanSTARRS lightcurves
import useful_funcs
import math
from math import log10, radians, pi,cos,sin
import plot

class lc:

    def __init__(self):
        
        self.lc_dir = "lightcurves/"
        useful_funcs.create_dir(self.lc_dir)
        self.query      = query.query_DES() # to get DES lightcurves
        self.query_S82  = query.Stripe82() # get SDSS S82 lightcurves
        self.query_SDSS = query.query_SDSS() # get SDSS spectra
        self.query_ZTF  = query.query_ZTF() # get ZTF lightcurves
        self.save_dir = self.lc_dir+"DES/"
        useful_funcs.create_dir(self.save_dir)
        self.quasar_catalog_dir = "catalog/"
        self.quasar_catalog = "DR14+DR7+OzDES+Milliq_S1S2.txt"
        self.lc_info_file = "lc_info.csv"
        self.log = "log_lightcurves"
        self.band_list = ["g","r","i","z"]
        self.save_line = ""

    ###########
    # General #
    ###########

    def set_output_dir(self,output_dir):

        # set the output dir (e.g. DES/, SDSS/)

        self.save_dir = output_dir
        useful_funcs.create_dir(self.save_dir)

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

    def fill_coadd_mag_with_median_mag(self,quasar_info,total_quasars):

        # get median magnitude for the quasar with wierd/without coadd mag 

        band = "i"
        if not ((quasar_info["mag_psf_%s" % band] > 0) & \
                (quasar_info["mag_psf_%s" % band] < 40 )):
            mag = np.median(total_quasars[total_quasars["band"]==\
                            band]["mag_psf"])
            quasar_info["mag_psf_%s" % band] = mag
        return quasar_info

    def write_header(self,quasar_catalog):

        if not os.path.isfile(self.lc_dir+self.lc_info_file):
            f = open(self.lc_dir+self.lc_info_file,"w")
            headers = ",".join(quasar_catalog.dtype.names)
            f.write("name,"+headers+",N_DES_g,N_DES_r,N_DES_i,N_DES_z")
            f.write(",N_SDSS_g,N_SDSS_r,N_SDSS_i,N_SDSS_z,N_PS_g,N_PS_r")
            f.write(",N_PS_i,N_PS_z,N_ZTF_g,N_ZTF_r\n")
            f.close()

    def save_information(self):

        f = open(self.lc_dir+self.lc_info_file,"a")
        f.write(self.save_line+"\n")
        f.close()

    def read_lc_info(self):

        lc_info = np.genfromtxt(self.lc_dir+self.lc_info_file,dtype=None,\
                                delimiter=",",names=True)
        return lc_info

    def get_unprocessed_quasars(self,quasar_catalog):    

        ra_total = quasar_catalog["ra"]
        subset= self.read_lc_info()
        ra_sub = subset["ra"]
        mask = np.array([tot not in ra_sub for tot in ra_total ])
        quasar_catalog = quasar_catalog[mask]

        return quasar_catalog

    ############
    # DES part #
    ############

    def generate_finalcut_lightcurve(self,quasar):

        # obtain the information at given position of the quasar.

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


    def convert_flux_to_mag(self,total_quasars):

        # convert flux to magnitude

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

    def save_DES_lightcurves(self,total_quasars):

        # save DES lightcurve in each band and record number of data points

        for band in self.band_list:
            data_in_band = total_quasars[total_quasars["band"]==band]
            if len(data_in_band) == 0 :
                print("No data found in "+band+" !")
                self.save_line = self.save_line+",0"
            else:
                self._process_and_save(data_in_band["mjd_obs"],\
                                       data_in_band["mag_psf"],\
                                       data_in_band["mag_err_psf"],band,\
                                       outlier=True,bin_data=True)

    #############
    # SDSS part #
    #############

    def generate_SDSS_lightcurves(self,quasar):

        # make the SDSS light curves from querying the database

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

    def save_SDSS_lightcurves(self,SDSS_quasars,mag_diff=None):

        # save the SDSS light curves into csv file

        for band in self.band_list:
            data_in_band = SDSS_quasars[SDSS_quasars["band"]==band]
            if len(data_in_band) == 0 :
                print "No data found in SDSS "+band+" !"
                self.save_line = self.save_line+",0"
            else:
                self._process_and_save(data_in_band["mjd_obs"],\
                                       data_in_band["mag_psf"],\
                                       data_in_band["mag_err_psf"],band,\
                                       outlier=True,bin_data=True,\
                                       mag_diff=mag_diff)

    #################
    # PanSTARRS DR2 #
    #################

    def generate_and_save_PS_lightcurves(self,quasar):

        # make and save PanSTARRS light curves

        dcolumns = ("""objID,detectID,filterID,obsTime,ra,dec,psfFlux,psfFluxErr,psfMajorFWHM,psfMinorFWHM,psfQfPerfect,apFlux,apFluxErr,infoFlag,infoFlag2,infoFlag3""").split(',')
        attempts = 0
        while attempts < 3:
            try:
                dresults = query_PS.ps1cone(quasar["ra"],quasar["dec"],2./3600.,\
                           table='detection',release='dr2',columns=dcolumns)
                break
            except requests.exceptions.HTTPError:
                attempts += 1
            time.sleep(60)
        if len(dresults) == 0: 
            print ("0 data points in Panstarrs ALL bands")
            self.save_line = self.save_line+",0,0,0,0"
            return 0
        dtab = query_PS.addfilter(ascii.read(dresults))
        if len(dtab) > 1:
            dtab.sort('obsTime')

        # save light curves

        for band in self.band_list:
            dtab_band = dtab[dtab['filter']==band]
            time = dtab_band['obsTime']
            mag = -2.5*np.log10(dtab_band['psfFlux']) + 8.90
            mag_error = dtab_band['psfFluxErr']/dtab_band['psfFlux']*1.09
            signal = 10**((22.5-mag)/2.5)
            error = mag_error*signal/1.09
            if len(dtab_band) == 0:
                print ("0 data points in Panstarrs "+band+" band !!")
                self.save_line = self.save_line+",0"
            else:
                self._process_and_save(time,mag,mag_error,band)

        return 0 


    ############
    # ZTF part #
    ############

    def generate_and_save_ZTF_lightcurves(self,quasar):

        data = self.query_ZTF.get_lightcurve(quasar["ra"],\
               quasar["dec"],2./3600.)
        if data is None or data.size == 0 :
            print ("0 data points in ZTF ALL bands")
            self.save_line = self.save_line+",0,0"
            return 0
        if data.size > 1:
            data.sort(order='mjd')
        for band in self.query_ZTF.band_list:
            data_band = data[data['filtercode'] == \
                        self.query_ZTF.band_name[band]]
            time = data_band["mjd"]
            mag = data_band["mag"]
            mag_error = data_band["magerr"]
            if (data_band) == 0:
                print("0 data points in ZTF "+band+" band !!")
                self.save_line = self.save_line+",0"
            else:
                self._process_and_save(time,mag,mag_error,band)
            print self.query_ZTF.band_list
        return 0

    ###########################
    # Pre-process lightcurves #
    ###########################

    def _process_and_save(self,time,signal,error,band,mag_diff=None,\
                          outlier=False,bin_data=False,record=True):
        # process the lightcurve (outlier rejection or binning) and save to 
        # the lightcurve directory

        # remove outlier 
        if outlier:
            mask = self._get_mask_sigma_clip_moving_avg(signal)
            time   = time[mask]
            signal = signal[mask]
            error  = error[mask]

        # bin the data at the same date
        if bin_data:
            time,signal,error = self._bin_data(time,signal,error)

        # save lightcurves
        
        np.savetxt("%s/%s.csv" % (self.save_dir , band),\
                   np.array([time,signal,error]).T,\
                   fmt="%f,%f,%f",comments="",\
                   header="mjd_obs,mag_psf,mag_err_psf")
        if record :
            self.save_line = self.save_line+","+str(len(signal))

        if mag_diff is not None: # For SDSS only
            corr_dir = self.save_dir.replace("SDSS", "SDSS_corr")
            useful_funcs.create_dir(corr_dir)
            signal = signal+mag_diff[band]
            np.savetxt("%s/%s.csv" % (corr_dir , band),\
                       np.array([time,signal,error]).T,\
                       fmt="%f,%f,%f",comments="",\
                       header="mjd_obs,mag_psf,mag_err_psf")


    def _get_mask_sigma_clip_moving_avg(self,signal,sigma_level=3):

        # return mask for sigma clip using moving median
        mask = np.ones( len(signal) , dtype=bool)
        outliers = True
        while outliers:
            sigma = np.std(signal[mask])
            signal_smooth = medfilt(signal[mask],kernel_size=5)
            mask_thistime = (abs(signal[mask] - signal_smooth) < \
                             sigma*sigma_level)
            if np.sum(~mask_thistime) == 0: outliers = False
            else: mask[mask] = mask_thistime
        return mask

    def _bin_data(self,time,signal,error):

        # bin the data points within the same day
        time2   = []
        signal2 = []
        error2  = []
        count   = 0

        while(count < len(time)):
            idx = ( np.floor(time) == np.floor(time[count]) )
            signal_temp = signal[idx]
            error_temp  = error[idx]
            nn          = len(signal_temp)

            signal_temp, error_temp = self.__mag2flux(signal_temp, error_temp)
            signal_temp, error_temp = self.__weighted_mean(signal_temp, error_temp)
            signal_temp, error_temp = self.__flux2mag(signal_temp, error_temp)

            time2.append( np.average(time[idx]) )
            signal2.append( signal_temp )
            error2.append( error_temp )

            count += nn

        time   = np.asarray(time2)
        signal = np.asarray(signal2)
        error  = np.asarray(error2)

        return time,signal,error

    def __mag2flux(self, signal, error):
        flux = 10.**(-1.*signal/2.5)
        return 10.**(-1.*signal/2.5), np.abs( -flux*error*np.log(10.)/2.5 )


    def __flux2mag(self, signal, error):
        return -2.5*np.log10(signal), np.abs( -2.5* error/signal/np.log(10.))


    def __weighted_mean(self, signal, error):
        signal_mean = np.sum(signal/error**2.) / np.sum(1./error**2.)
        error_mean  = np.sqrt( np.sum(error**2.) ) / np.sqrt( np.float(len(signal)) )
        return signal_mean, error_mean


    ###################
    # lightcurve stat #
    ###################

    def get_clean_sample(self,quasars_info,N_DES=50,N_SDSS=30,N_band=2):

        print ("Parent_sample: "+str(len(quasars_info)))
        N = 0
        clean_sample = np.array([],dtype=quasars_info.dtype)
        for row in quasars_info:
            N_pass = 0
            for band in self.band_list:
                if( row["N_DES_"+band]>N_DES) and (row["N_SDSS_"+band]>N_SDSS) \
                    and (row["spread_model_i"] < 0.005):
                    N_pass = N_pass +1
            if N_pass > N_band-1:
                clean_sample = np.append(clean_sample,row)
        print ("Number of Quasar with Enough epochs: "+str(len(clean_sample)))
        print ("z:"+str(np.median(clean_sample["z"])))
        print ("mag_psf_i:"+str(np.median(clean_sample["mag_psf_i"])))
        print ("N_DES_i:"+str(np.median(clean_sample["N_DES_i"])))
        print ("N_SDSS_i:"+str(np.median(clean_sample["N_SDSS_i"])))

        np.savetxt(self.lc_dir+"lc_clean.csv",clean_sample,\
                   fmt="%s,"*(len(clean_sample.dtype.names)-1)+"%s",\
                   header=",".join(clean_sample.dtype.names),comments="")

        return clean_sample

############################################################################
    ##################
    # Redudant Funcs #
    ##################

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

############################################################################

class spectra:

    def __init__(self):

        self.query_SDSS = query.query_SDSS()
        self.dir_filter = "filters/"  # create filter info dir
        useful_funcs.create_dir(self.dir_filter)
        self.dir_spec = "spectra/" # create dir for spectra
        useful_funcs.create_dir(self.dir_spec)
        self.band_list = ["g","r","i","z"]
        self.filter_info = {"SDSS":{},"DES":{},"LCO":{}}
        self.get_DES_SDSS_bandpass()
        self.get_LCO_bandpass()


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
            self.filter_info["SDSS"].update({band:data_SDSS[:,[0,1]]})

        # load DES bandpass info
        band_num = {"g":1,"r":2,"i":3,"z":4}
        data_DES = np.genfromtxt(self.dir_filter+"DES.dat")
        for band in self.band_list:
            self.filter_info["DES"].update({band:data_DES[:,[0,band_num[band]]]})

    def get_LCO_bandpass(self):

        # dowload SDSS bandpass
        url_filter_LCO = "https://lco.global/documents/%s/sdss.%sp.txt"
        band_number = {"g":10,"r":11,"i":12}
        for band in self.band_list[0:3]:
            useful_funcs.download_file(url_filter_LCO % (band_number[band],band),
                                       self.dir_filter+"LCO_"+band+".dat")
        for band in self.band_list[0:3]:
            data_LCO = np.genfromtxt(self.dir_filter+"LCO_"+band+".dat",\
                                     skip_header=1)
            data_LCO[:,0] = data_LCO[:,0]*10
            self.filter_info["LCO"].update({band:data_LCO})

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

    def mag_LCO_to_DES(self,name):

        mag_diff = {}
        if not os.path.isfile(self.dir_spec+name+".fits"):
            for band in self.band_list[0:3] : mag_diff.update({band:0.0})
            return mag_diff

        hdul = fits.open(self.dir_spec+name+".fits")
        data =  hdul[1].data
        data_clean = data[(data["ivar"]>0) & (data["and_mask"] == 0)]

        lambda_q = 10**data_clean["loglam"]
        flux_q =  10**(-17)*data_clean["flux"]
        for band in self.band_list[0:3] :
            trans_LCO,trans_DES = self.filter_info["LCO"][band],\
                                   self.filter_info["DES"][band]
            mag_LCO =  self.spectrum_to_mag(lambda_q,flux_q,trans_LCO[:,0],\
                                             trans_LCO[:,1])
            mag_DES =  self.spectrum_to_mag(lambda_q,flux_q,trans_DES[:,0],\
                                            trans_DES[:,1])
            #print ("mag_SDSS: "+str(mag_SDSS)+" ; mag_DES: "+str(mag_DES))
            diff = mag_DES-mag_LCO
            print (band+": m_DES - m_LCO =" + str(diff))
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



