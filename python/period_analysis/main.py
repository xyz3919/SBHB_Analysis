import numpy as np
import pandas as pd
import os
import sys
from quasar_drw import quasar_drw as qso_drw
from plot import plot
import useful_funcs
from scipy import stats
from scipy.optimize import curve_fit

class analysis:

    def __init__(self):

        self.lc_dir = "lightcurves/"
        self.lc_info_file = "lc_clean.csv"
        self.output_dir = "analysis/"
        useful_funcs.create_dir(self.output_dir)
        self.random_state = np.random.RandomState(0)
        self.band_list = ["g","r","i","z"]
        self.surveys = ["DES","SDSS_corr"]
        self.lc_names = ["mjd_obs","mag_psf","mag_err_psf"]
        self.period_lowerlim = 500

    def check_lightcurves_exist(self,name,band):

        exist = True
        for survey in self.surveys:
            file_path = self.lc_dir+survey+"/"+name+"/"+band+".csv"
            exist = exist and os.path.exists(file_path)

        return exist

    def read_quasar_catalog(self):

        df = pd.read_csv(self.lc_dir+self.lc_info_file)
        return df

    def read_lightcurve(self,name,survey,band):

        file_path = self.lc_dir+survey+"/"+name+"/"+band+".csv"
        data =  pd.read_csv(file_path,comment="#")
        return data
    def error_boostraping(self,lc,periodogram,band):

        number_times = 100
        period_obs, signal_obs = lc.periodogram(lc.time,lc.signal)
        signal_boost = self.random_state.normal(lc.signal, lc.error, \
                       [number_times,len(lc.signal)])
        signal_total = []
        for signal in signal_boost:
            period_boost, psd_boost = lc.periodogram(lc.time,signal)
            signal_total.append(psd_boost)

        signal_total = np.array(signal_total, dtype=np.float64)

        error_list = []
        for i in range(len(signal_total[0])):
            error = np.std(signal_total[:,i])
            error_list.append(error)
        error_list = np.array(error_list, dtype=np.float64)
        periodogram.plot_boost_periodogram(period_obs, signal_obs, error ,band)

    def check_period_max_amp(self,_freq, psd):
        # only search for peak amp between 500 days and 3 cycles.
        try:
            period = _freq[_freq<max(_freq)/3.][_freq>self.period_lowerlim]
            amp = psd[_freq<max(_freq)/3.][_freq>self.period_lowerlim]
            period_max = period[np.where(amp==np.max(amp))]
            return period_max

        except: return None

    def fitting(self,time,signal,error,period):

        def sin_func(x,amplitude,ref_day,median):
            return amplitude*np.sin(2*np.pi*x/period+ref_day)+median
        p0 = [1,50000,20]
        popt, pcov = curve_fit(sin_func, time, signal,p0=p0)

        xn = np.linspace(np.min(time)-100,np.max(time)+100,10000)
        yn = sin_func(xn,*popt)
        return xn,yn

    def clean_parameters_list(self,parameters_list):

        num_parameters = len(parameters_list[0])
        for i in range(0,num_parameters-1):
            array = parameters_list[:,i]
            upper = np.percentile(array,90)
            lower = np.percentile(array,10)
            parameters_list = parameters_list[(array<upper) & (array>lower)]
        return parameters_list

    def save_period_amp(self,_freq, psd,filename):

        save_data = zip(_freq,psd)
        np.savetxt(filename,save_data,delimiter=",",header="period,amplitude")

    def save_period_confidence_level(self,_freq,boundary_all,filename):

        period = np.array([_freq]).T
        psd_array = np.array([boundary_all]).T
        save_data = np.concatenate((period, psd_array), axis=1)
        np.savetxt(filename,save_data,delimiter=",",header="period,confidence_level")

    def tailored_simulation(self,lc,band,z,name,periodogram,lightcurve):

        psd_mock_all = []
        nwalkers = 10
        burnin = 10
        Nsteps = 20
        draw_times = 100
        parameters_list =  lc.fit_drw_emcee(nwalkers=nwalkers, burnin=burnin,\
                           Nstep=Nsteps,random_state=self.random_state)
        parameters_list_good = self.clean_parameters_list(parameters_list)
        for i in range(draw_times):
            tau,c,b = np.exp(parameters_list_good[self.random_state.randint(\
                      len(parameters_list_good))])
            mock_time,mock_signal = lc.generate_mock_lightcurve(tau,b,c,lc.time,\
                                    z,random_state=self.random_state)
            #lightcurve.plot_mock_curve(mock_time,mock_signal_correct,band)
            period_mock, psd_mock = lc.periodogram(mock_time,mock_signal)
            psd_mock_all.append(psd_mock)
            periodogram.plot_mock_periodogram(period_mock, psd_mock,band)

        psd_at_each_period = zip(*psd_mock_all)
        confidence_levels = []
        period_obs, psd_obs = lc.ls_astroML()
        for i in range(len(period_obs)):
            confidence_level_at_each_period = float(stats.percentileofscore(\
                                              psd_at_each_period[i],psd_obs[i]))
            confidence_levels.append(100.-confidence_level_at_each_period)
        self.save_period_confidence_level(period_obs,confidence_levels,\
                                   self.output_dir+name+"/confidence_"+\
                                   band+".csv")
        periodogram.plot_confidence_level(period_mock, psd_mock_all, band)

    def analyze_lightcurve(self,quasar):

        periodogram = plot(2,2)
        lightcurve = plot(4,1,figsize=(8,8),sharex=True)

        name = quasar["name"]
        print ("-- Analyzing quasar "+name+ "--")
        useful_funcs.create_dir(self.output_dir+name)

        for band in self.band_list:
            if not self.check_lightcurves_exist(name,band):
                print ("SDSS and DES lightcurves not found in "+band+" band !!")
            else:
                time_all,signal_all,error_all = [],[],[]
                # read main survey 
                survey = self.surveys[0]
                data = self.read_lightcurve(name,survey,band)
                lc = qso_drw(data["mjd_obs"],data["mag_psf"],\
                             data["mag_err_psf"], quasar["z"],\
                             preprocess=True)
                time, signal, error = lc.get_lc()
                lightcurve.plot_light_curve(time,signal,error,survey,band)

                # add other surveys
                for survey in self.surveys[1:]:
                    data = self.read_lightcurve(name,survey,band)
                    lc2 = qso_drw(data["mjd_obs"],data["mag_psf"],\
                                 data["mag_err_psf"], quasar["z"],\
                                 preprocess=True)
                    time, signal, error = lc2.get_lc()
                    lightcurve.plot_light_curve(time,signal,error,survey,band)
                    lc.add_lc(time, signal, error,preprocess=False)

                # periodogram
                if not len(lc.time)>3:
                    print ("Not enough data points !!")
                else:
                    period, psd = lc.ls_astroML()
                    periodogram.plot_periodogram(period, psd,band)
                    self.error_boostraping(lc,periodogram,band)

                    period_max = self.check_period_max_amp(period, psd)
                    if period_max is  None:
                        print ("Can not find the peak value !!")
                    else:
                        xn,yn = self.fitting(time,signal,error,period_max)
                        lightcurve.plot_fit_curve(xn,yn,band)

                    self.tailored_simulation(lc,band,quasar["z"],name,\
                                             periodogram,lightcurve)
                    self.save_period_amp(period, psd,self.output_dir+name+\
                                       "/periodogram_"+band+".csv")
        lightcurve.savefig(self.output_dir+name,"/lightcurve.png",name)
        periodogram.savefig(self.output_dir+name,"/periodogram.png",name)





                    

                    



