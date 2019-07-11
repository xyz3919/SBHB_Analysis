import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('agg')
from javelin.zylc import get_data
from javelin.lcmodel import Cont_Model
from gatspy import datasets, periodic
import carmcmc as cm
from quasar_drw import quasar_drw as qso_drw
from quasar_drw import lnlike
from plot import plot
import useful_funcs

class analysis:

    def __init__(self):

        self.lc_dir = "lightcurves/"
        self.lc_info_file = "lc_clean.csv"
        #self.lc_info_file = "lc_remain.csv"
        self.output_dir = "analysis/"
        useful_funcs.create_dir(self.output_dir)
        self.stat_dir = "statistics/"
        useful_funcs.create_dir(self.stat_dir)
        self.carma_dir = "carma/"
        self.random_state = np.random.RandomState(0)
        self.band_list = ["g", "r", "i", "z"]
        self.surveys = ["DES", "SDSS_corr"]
        self.surveys_add = ["ZTF","PS"]
        self.lc_names = ["mjd_obs", "mag_psf", "mag_err_psf"]
        self.period_lowerlim = 500
        self.test = False

    def check_lightcurves_exist(self, name, band):

        exist = True
        for survey in self.surveys:
            file_path = self.lc_dir+survey+"/"+name+"/"+band+".csv"
            exist = exist and os.path.exists(file_path)

        return exist

    def make_combined_lc_for_javelin(self, lc, name, band):

        zipdata = zip(lc.time,lc.signal,lc.error)
        filename = self.lc_dir+"combined/"+name+"/"+band+".dat"
        np.savetxt(filename,zipdata,delimiter=" ",comments="",fmt='%f')

        return filename


    def read_quasar_catalog(self):

        #lc_info = pd.read_csv(self.lc_dir+self.lc_info_file)
        lc_info = np.genfromtxt(self.lc_dir+self.lc_info_file, delimiter=",",\
                  skip_header=1,\
                  dtype=[("name", "|S20"), ("ra", float), ("dec", float),\
                         ("z", float), ("flag_0", "|S6"), ("flag_1", "|S6"),\
                         ("flag_2", "|S6"), ("flag_3", "|S15"),\
                         ("mag_i", float), ("spread_model_i", float),\
                         ("spread_model_err_i", float),\
                         ("N_DES_g", int), ("N_DES_r", int), ("N_DES_i", int),\
                         ("N_DES_z", int), ("N_SDSS_g", int), ("N_SDSS_r", int),\
                         ("N_SDSS_i", int), ("N_SDSS_z", int)])
        return lc_info

    def read_lightcurve(self, name, survey, band):

        file_path = self.lc_dir+survey+"/"+name+"/"+band+".csv"
        if not os.path.exists(file_path):
            return None
        else:
            data = pd.read_csv(file_path, comment="#")
            return data

    def error_boostraping(self, lc, band):

        if self.test is True: number_times = 100
        else: number_times = 1000
        period_obs, signal_obs = lc.periodogram(lc.time, lc.signal)
        signal_boost = self.random_state.normal(lc.signal, lc.error, \
                       [number_times, len(lc.signal)])
        signal_total = []
        for signal in signal_boost:
            period_boost, psd_boost = lc.periodogram(lc.time, signal)
            signal_total.append(psd_boost)

        #np.save("%s/psds_boost_%s.npy" % (save_dir,band),signal_total)

        signal_total = np.array(signal_total, dtype=np.float64)

        error_list = []
        for i in range(len(signal_total[0])):
            error = np.std(signal_total[:, i])
            error_list.append(error)
        error_list = np.array(error_list, dtype=np.float64)
        #periodogram.plot_boost_periodogram(period_obs, signal_obs, error, band)
        
        return error_list

    def _check_period_max_amp(self, period, psd):
        # only search for peak amp between 500 days and 3 cycles.
        interval = (period<max(period)/3.) & (period>self.period_lowerlim)
        period_within = period[interval]
        amp_within = psd[interval]
        amp_max = np.max(amp_within)
        period_max = period_within[np.where(amp_within==amp_max)]
        return period_max,amp_max


    def _fitting(self, time, signal, period):

        def sin_func(x, amplitude, ref_day, median):
            return amplitude*np.sin(2*np.pi*x/period+ref_day)+median
        p0 = [1, 50000, 20]
        popt, pcov = curve_fit(sin_func, time, signal, p0=p0)

        xn = np.linspace(np.min(time)-100, np.max(time)+100, 10000)
        yn = sin_func(xn, *popt)
        return xn, yn

    def _exp_sin(self,x,T,tau):
        return np.cos(2*np.pi/T*x)*np.exp(-x/tau)

    def _fitting_exp_sin(self, time, signal):

        p0 = [3 , 1 ]
        popt, pcov = curve_fit(self._exp_sin,time,signal, p0=p0, \
                               bounds=(0,[20,20]),maxfev=5000)
        # calculate goodness of fit
        residuals = signal- self._exp_sin(time, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((signal-np.mean(signal))**2)
        r_squared = 1 - (ss_res / ss_tot)

        perr = np.sqrt(np.diag(pcov))
        xn = np.linspace(np.min(time)-0.1, np.max(time)+0.1, 10000)
        yn = self._exp_sin(xn, *popt)
        return popt,perr,xn,yn,r_squared

    def clean_parameters_list(self, parameters_list):

        num_parameters = len(parameters_list[0])
        excluded = np.zeros(len(parameters_list[:,0]), dtype=bool)
        for i in range(0, num_parameters-1):
            array = parameters_list[:, i]
            upper = np.percentile(array, 95)
            lower = np.percentile(array, 5)
            excluded = (excluded) | ( ~((array<upper) & (array>lower) ) ) 
        parameters_list = parameters_list[~excluded]
        return parameters_list

    def save_lightcurve(self,time,signal,error,filename):

        save_data=zip(time,signal,error)
        np.save(filename,save_data)
        #np.savetxt(filename, save_data, delimiter=",", comments="",\
        #           fmt="%1.4E",header="time,signal,error")

    def save_period_amp(self, _freq, psd, error, filename):

        if error == None:
            save_data = zip(_freq, psd)
            np.save(filename,save_data)
            #np.savetxt(filename, save_data, delimiter=",", comments="",\
            #           fmt="%1.4E",header="period,amplitude")
        else:
            save_data = zip(_freq, psd, error)
            np.save(filename,save_data)
            #np.savetxt(filename, save_data, delimiter=",", comments="",\
            #           fmt="%1.4E",header="period,amplitude,error")

    def save_drw_parameters(self, tau, var, mu,filename):

        # save DRW parameters (tau, variance and mean) 

        save_data = zip(tau,var,mu)
        np.savetxt(filename, save_data, delimiter=",", comments="",\
                   header="tau,var,mu")

    def add_flux_fluxerr_into_data(self,data):

        """ add flux_psf and flux_err_psf into data """

        data["flux_psf"] = 10**((22.5-data["mag_psf"])/2.5)
        data["flux_err_psf"] = data["mag_err_psf"]*data["flux_psf"]/1.09

        return data

    def save_period_confidence_level(self, _freq, boundary_all, filename):

        period = np.array([_freq]).T
        psd_array = np.array([boundary_all]).T
        save_data = np.concatenate((period, psd_array), axis=1)
        np.savetxt(filename, save_data, delimiter=",", comments="",\
                   header="period,confidence_level")

    def tailored_simulation(self,lc,band,z,name):#,periodogram,lightcurve):

        """ running the DRW model fitted to each quasar """

        mock_dir = self.output_dir+name+"/mock"
        useful_funcs.create_dir(mock_dir)

        #psd_mock_all = []
        # parameters
        if self.test is True: 
            nwalkers,burnin,Nsteps,draw_times = 100,50,100,100
        else:
            nwalkers = 500
            burnin = 150
            Nsteps = 500
            draw_times = 5000

        # fitting DRW model
        samples =  lc.fit_drw_emcee(nwalkers=nwalkers, burnin=burnin,\
                           Nstep=Nsteps,random_state=self.random_state)

        # plot walkers
        walkers = plot(3,1,figsize=(10,7), sharex=True)
        walkers.plot_walkers(np.exp(samples))
        walkers.savefig(mock_dir,"/walker_"+band+".png",band+" band")

        parameters_list = samples[:, burnin:, :].reshape((-1, 3))

        # remove top 5% and bottome 5%
        parameters_list_good = self.clean_parameters_list(parameters_list)

        # plot posterior
        theta = [parameters_list_good[:,0],parameters_list_good[:,1],\
                 parameters_list_good[:,2]]
        likelihood = []
        for theta in parameters_list_good:
            likelihood.append(lnlike(theta,lc.time,lc.signal,lc.error,z))
        from plot import plot_posterior
        plot_posterior(np.exp(parameters_list_good),likelihood,\
                       band,mock_dir+"/post_%s.png" % band)

        # make the mock light curves from DRW parameters

        mock_lcs,mock_psds = [],[]
        for i in range(draw_times):
            tau,c,b = np.exp(parameters_list_good[self.random_state.randint(\
                    len(parameters_list_good))])
            mock_time,mock_signal = lc.generate_mock_lightcurve(tau,c,lc.time,\
                                    lc.signal,z,random_state=self.random_state)
            mock_lcs.append(mock_signal)
            #mock_lcs.append(np.array([mock_time,mock_signal]).T)
            #self.save_lightcurve(mock_time,mock_signal,lc.error,\
            #                     "%s/lc_%s.csv" % (mock_dir,i))
            #lightcurve.plot_mock_curve(mock_time,mock_signal_correct,band)

            period_mock, psd_mock = lc.periodogram(mock_time,mock_signal)
            mock_psds.append(psd_mock)
            #self.save_period_amp(period_mock,psd_mock,None,\
            #                     "%s/psd_%s.csv" % (mock_dir,i))
            #mock_psds.append(np.array([period_mock,psd_mock]).T)
            #periodogram.plot_mock_periodogram(period_mock, psd_mock,band)
        np.save("%s/lcs_%s.npy" % (mock_dir,band),mock_lcs)
        np.save("%s/psds_%s.npy" % (mock_dir,band),mock_psds)

        """
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
        """

    def tailored_simulation_javelin(self,lc,band,z,name,periodogram,lightcurve):

        useful_funcs.create_dir(self.lc_dir+"combined")
        useful_funcs.create_dir(self.lc_dir+"combined/"+name)

        filename = self.make_combined_lc_for_javelin(lc,name,band)
        javdata = get_data(filename,names=["Continuum"])
        cont = Cont_Model(javdata)
        cont.do_mcmc(nwalkers=500, nburn=100, nchain=500,fchain="test.dat")
        cont.show_hist(figout="javelin_%s" % band, figext=".png")

    def do_multi_band_periodogram(self,lightcurves_total,multi_periodogram):

        time,flux,flux_err,band = lightcurves_total
        time = time.astype(float)
        flux = flux.astype(float)
        flux_err = flux_err.astype(float)
        T_max = float(np.max(time)-np.min(time))
        T_min = float(T_max/len(flux))
        model = periodic.LombScargleMultiband(fit_period=True)
        model.optimizer.period_range=(T_min,T_max)
        model.fit(time,flux,flux_err,band)
        periods = np.linspace(T_min,T_max,10*len(flux))
        P_multi = model.periodogram(periods)
        multi_periodogram.plot_multi_periodogram(periods,P_multi,"total")

    def do_carma_process(self,lc,band,name,periodogram,lightcurve):

        """ Using carma-pack to fit DRW model"""

        p = 1
        q = 0
        sample_number = 100000
        model = cm.CarmaModel(lc.time, lc.signal, lc.error, p=p, q=q)
        sample = model.run_mcmc(sample_number)
        tau = 1./np.exp(sample.get_samples("log_omega"))
        var = sample.get_samples("var")
        mu = sample.get_samples("mu")
        array = np.array([tau,var,mu]).T[0]
        from plot import plot_posterior_carma
        plot_posterior_carma(array,band,\
                             self.output_dir+name+"/post_carma_%s.png" % band)
        self.save_drw_parameters(tau,var,mu,self.output_dir+name+\
                                 "/post_carma_%s.txt" % band)

        npaths = 5000
        psd_mock_all = []
        for i in range(npaths):
            time_cont = np.linspace(min(lc.time),max(lc.time),\
                                    int((max(lc.time)-min(lc.time))/2.))
            time_cont_sim = time_cont+max(lc.time)
            signal_sim_cont = sample.simulate(time_cont_sim, bestfit='random')
            signal_sim = np.interp(lc.time, time_cont, signal_sim_cont)
            period_mock,psd_mock = lc.periodogram(lc.time,signal_sim)
            psd_mock_all.append(psd_mock)
            periodogram.plot_mock_periodogram(period_mock,psd_mock,band)
            #lightcurve.plot_mock_curve(time_cont,signal_sim_cont,band)
        periodogram.plot_confidence_level(period_mock,psd_mock_all,band)
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

    def _calculate_confidence_level(self,period,psd,psds_mock,band,ACF=False):

        psd_mock_at_each_period = zip(*psds_mock)
        confidence_levels = []
        for i in range(len(period)):
            confidence_level = float(stats.percentileofscore(\
                                     psd_mock_at_each_period[i],psd[i]))
            confidence_levels.append(100.-confidence_level)
        if ACF: append = "ACF" 
        else: append = ""
        self.save_period_confidence_level(period,confidence_levels,\
                                   self.output_dir+self.name+"/confidence_"+\
                                   append+"_"+band+".csv")
        return confidence_levels

    def _load_confidence_level(self,name,band,ACF=True):

        if ACF: append = "ACF"
        else: append = ""
        data = np.genfromtxt(self.output_dir+name+"/confidence_"+append+\
                             "_"+band+".csv",\
                             names=True,delimiter=",")
        return data

    def _get_fit_curve(self,time,signal,name,band):

        """ make fitted sin curve  """
        confidence_level = self._load_confidence_level(name,band)
        period_max,amp_max = self._check_period_max_amp(\
                             confidence_level["period"],\
                             100-confidence_level["confidence_level"])
        if len(period_max) > 1 : period_max = np.mean(period_max)
        xn, yn = self._fitting(time, signal, period_max)

        return xn,yn


    def plot_periodogram_and_lightcurve(self,quasar):

        """ making the periodogram and light curve plot """
        name = quasar["name"]
        self.name = name
        print ("-- Making periodogram plot %s -- " % name)

        real_dir = self.output_dir+name+"/real"
        mock_dir = self.output_dir+name+"/mock"

        periodogram = plot(2,2)
        for band in self.band_list:
            if os.path.exists("%s/psd_%s.npy" % (real_dir,band)):
                period,psd,psd_error = np.load("%s/psd_%s.npy" % \
                                       (real_dir,band))
                psd_mock = np.load("%s/psds_%s.npy" % (mock_dir,band))
                significance_level = self._calculate_confidence_level(\
                                     period,psd,psd_mock,band)
                periodogram.plot_periodogram(period, psd,band)
                periodogram.plot_boost_periodogram(period, psd, psd_error, band)
                periodogram.plot_confidence_level(period, psd_mock, band)
                periodogram.plot_mock_periodogram(period, psd_mock, band)
                periodogram.plot_peak_period(period, significance_level,band)
        periodogram.savefig(self.output_dir+name,"/periodogram.png",name)

        print ("-- Making light curve plot %s --" % name)

        lightcurve = plot(4,1,figsize=(8,8),sharex=True)
        for band in self.band_list:
            for survey in self.surveys+self.surveys_add:
                data = self.read_lightcurve(name,survey,band)
                if data is not None:
                    time, signal, error = data["mjd_obs"],data["mag_psf"],\
                                          data["mag_err_psf"]
                    if survey in self.surveys_add: adjust_lim = False
                    else: adjust_lim = True
                    lightcurve.plot_light_curve(time,signal,error,survey,band,\
                                                adjust_lim=adjust_lim)

            if os.path.exists("%s/lc_%s.npy" % (real_dir,band)):
                time,signal,error = np.load("%s/lc_%s.npy" % \
                                    (real_dir,band))
                signal_mock = np.load("%s/lcs_%s.npy" % (mock_dir,band))
                xn, yn = self._get_fit_curve(time,signal,name,band)
                lightcurve.plot_fit_curve(xn,yn,band)
        lightcurve.savefig(self.output_dir+name,"/lightcurve.png",name)


    def analyze_lightcurve(self,quasar):

        """ the main analyzing function """

        name = quasar["name"]
        print ("-- Analyzing quasar %s --" % name)

        useful_funcs.create_dir(self.output_dir+name)
        real_dir = self.output_dir+name+"/real"
        useful_funcs.create_dir(real_dir)

        for band in self.band_list:
            # avoid missing band
            if not self.check_lightcurves_exist(name,band):
                print ("SDSS and DES lightcurves not found in "+band+" band !!")
            else:
                # read DES survey
                survey = self.surveys[0]
                data = self.read_lightcurve(name,survey,band)
                #data = self.add_flux_fluxerr_into_data(data)
                time, signal, error = data["mjd_obs"],data["mag_psf"],\
                                      data["mag_err_psf"]
                lc = qso_drw(time, signal, error, quasar["z"], preprocess=False)
                # add other surveys
                for survey in self.surveys[1:]:
                    data = self.read_lightcurve(name,survey,band)
                    #data = self.add_flux_fluxerr_into_data(data)
                    if len(data) > 0 :
                        time,signal,error = data["mjd_obs"],data["mag_psf"],\
                                            data["mag_err_psf"]
                        lc.add_lc(time, signal, error, preprocess=False)
                np.save("%s/lc_%s.npy" % (real_dir,band),\
                        [lc.time,lc.signal,lc.error])
                
                # periodogram
                if not len(lc.time)>3:
                    print ("Not enough data points !!")
                else:
                    # run Lomb-Scargle for target
                    period, psd = lc.ls_astroML()
                    # save error estimated from boostrping
                    psd_error = self.error_boostraping(lc,band)

                    np.save("%s/psd_%s.npy" % (real_dir,band),\
                            [period, psd,psd_error])
                    
                    # tailored simulation
                    self.tailored_simulation(lc,band,quasar["z"],name)

    def analyze_lightcurve_orig(self,quasar):

        """ the main analyzing function """

        #periodogram = plot(2,2)
        #multi_periodogram = plot(1,1)
        #periodogram_carma = plot(2,2)
        #lightcurve = plot(4,1,figsize=(8,8),sharex=True)

        name = quasar["name"]
        print ("-- Analyzing quasar "+name+ "--")
        useful_funcs.create_dir(self.output_dir+name)
        #useful_funcs.create_dir(self.output_dir+name+"/"+self.carma_dir)
        lightcurves_total = []

        for band in self.band_list:
            if not self.check_lightcurves_exist(name,band):
                print ("SDSS and DES lightcurves not found in "+band+" band !!")
            else:
                time_all,signal_all,error_all = [],[],[]
                # read main survey 
                survey = self.surveys[0]
                data = self.read_lightcurve(name,survey,band)
                data = self.add_flux_fluxerr_into_data(data)
                lc = qso_drw(data["mjd_obs"],data["mag_psf"],\
                             data["mag_err_psf"], quasar["z"],\
                             preprocess=False)
                time, signal, error = lc.get_lc()

                #lightcurve.plot_light_curve(time,signal,error,survey,band)

                # add other surveys
                for survey in self.surveys[1:]:
                    data = self.read_lightcurve(name,survey,band)
                    data = self.add_flux_fluxerr_into_data(data)
                    lc2 = qso_drw(data["mjd_obs"],data["mag_psf"],\
                                  data["mag_err_psf"], quasar["z"],\
                                 preprocess=False)
                    time, signal, error = lc2.get_lc()
                    if len(signal) > 0 :
                        lightcurve.plot_light_curve(time,signal,error,\
                                                    survey,band)
                        lc.add_lc(time, signal, error,preprocess=False)

                # periodogram
                if not len(lc.time)>3:
                    print ("Not enough data points !!")
                else:
                    period, psd = lc.ls_astroML()
                    periodogram.plot_periodogram(period, psd,band)
                    multi_periodogram.plot_multi_periodogram(period,psd,band)
                    #periodogram_carma.plot_periodogram(period,psd,band)
                    error_list = self.error_boostraping(lc,periodogram,band)
                    #error_list = self.error_boostraping(lc,periodogram_carma,band)

                    period_max,amp_max = self._check_period_max_amp(period, psd)
                    if period_max is  None:
                        print ("Can not find the peak value !!")
                    else:
                        periodogram.plot_peak_period(period_max,band)
                        xn,yn = self._fitting(lc.time,lc.signal,lc.error,period_max)
                        lightcurve.plot_fit_curve(xn,yn,band)
                    lightcurves_total.append(np.array([lc.time,lc.signal,\
                                             lc.error,[band]*len(lc.time)])) 
                    self.tailored_simulation(lc,band,quasar["z"],name,\
                                             periodogram,lightcurve)
                    #self.do_carma_process(lc,band,name,periodogram_carma,lightcurve)
                    self.save_period_amp(period, psd,error_list,\
                                         self.output_dir+name+\
                                         "/periodogram_"+band+".csv")
        lightcurves_total = np.concatenate(lightcurves_total, axis=1)
        np.savetxt(self.output_dir+name+"/lightcurve_total.csv",\
                   lightcurves_total, delimiter=",", \
                   comments="",header="time,signal,error,band")
        self.do_multi_band_periodogram(lightcurves_total,\
                                       multi_periodogram)

        lightcurve.savefig(self.output_dir+name,"/lightcurve.png",name)
        periodogram.savefig(self.output_dir+name,"/periodogram.png",name)
        multi_periodogram.savefig(self.output_dir+name,"/periodogram_multi.png",name)
        #periodogram_carma.savefig(self.output_dir+name,"/periodogram_carma.png",name)

    ########################
    # Selecting candidates #
    ########################

    def record_confidence_peak(self):

        lower_period = self.period_lowerlim
        f_confidence = open(self.stat_dir+"candidates_confidence.txt","w")
        f_confidence.write("name,peak_g,peak_r,peak_i,peak_z\n")
        quasar_catalog = self.read_quasar_catalog()
        for quasar in quasar_catalog:
            name = quasar["name"]
            f_confidence.write(name)
            for band in self.band_list:
                if os.path.exists(self.output_dir+name+"/confidence_"+band+".csv"):
                    periodogram = pd.read_csv(self.output_dir+name+"/confidence_"+band+".csv",\
                                              names=["period","confidence"],skiprows=1)
                    upper_period =max(periodogram["period"])/3
                    periodogram1 = periodogram[periodogram["period"]<upper_period]
                    periodogram2= periodogram1[periodogram1["period"]>lower_period]
                    if len(periodogram2) != 0 :
                        peak_confidence = min(periodogram2["confidence"])
                        f_confidence.write(","+str(peak_confidence))
                    else: f_confidence.write(",-99")
                else: f_confidence.write(",999")
            f_confidence.write("\n")
        f_confidence.close()

    def find_strong_candidates(self):

        """ Compare the power with mock light curves """
    #    cand_strong = open("statistics/strong_candidates.txt","w")
        cand = pd.read_csv(self.stat_dir+"candidates_confidence.txt")
    #    for band in bands:
    #        cand = cand[cand["peak_"+band]<1]
        strong_cand = pd.DataFrame(columns=["name","peak_g","peak_r","peak_i","peak_z"])
        cand["yes"] = 0
        for index, row in cand.iterrows():
            yes = 0
            for band in self.band_list:
                if row["peak_"+band]<0.3:
                    yes= yes+1
            cand.at[index,"yes"] = yes
        strong_cand = cand[cand["yes"]>1]
        strong_cand.to_csv(self.stat_dir+"strong_candidates.csv",index=False)

    def generate_drw_parameters_plot(self):

        """ Plot the DRW parameters for the whole sample"""
        for band in self.band_list:
            print band
            quasar_catalog = self.read_quasar_catalog()
            tau_total = []
            var_total = []
            for quasar in quasar_catalog:
                name = quasar["name"]
                z = quasar["z"]
                print name
                if os.path.exists(self.output_dir+name+"/post_carma_%s.txt" %\
                                 band):
                    parameters = pd.read_csv(self.output_dir+name+\
                                 "/post_carma_%s.txt" % band,\
                                 names = ["tau","var","mu"],skiprows=1)
                    median = np.median(parameters,axis=0)
                    median[0] = median[0]/(1+z)
                    median = np.log10(median)
                    tau_total.append(median[0])
                    var_total.append(median[1])

            from plot import plot_drw_parameters
            plot_drw_parameters(tau_total,var_total,band,\
                                self.stat_dir+"drw_paramets_%s.png" % band)
         

    def make_pdf_for_strong_candidates(self):

        """ making the latex file and figure dir for stong cadidates """

        quasars = pd.read_csv(self.stat_dir+"strong_candidates.csv")
        header = "\documentclass[12pt]{article}\n"+\
                 "\\usepackage{graphicx}\n"+\
                 "\\usepackage[export]{adjustbox}\n"+\
                 "\\usepackage[margin=0.5in]{geometry}\n"+\
                 "\\begin{document}\n"
        section = "\\section*{%s}\n"
        body1 = "\\begin{figure}[!bp]\n"+\
                "\\begin{minipage}[b]{0.48\\textwidth}\n"+\
                "\\includegraphics[width=\\textwidth]"+\
                "{%s/periodogram.png}\n"+\
                "\\end{minipage}\n"+\
                "\\hfill\n"+\
                "\\begin{minipage}[b]{0.48\\textwidth}\n"+\
                "\\includegraphics[width=\\textwidth,right]"+\
                "{%s/lightcurve.png}\n"+\
                "\\end{minipage}\n"+\
                "\\includegraphics[width=0.8\\textwidth]"+\
                "{%s/ACF.png}\n"+\
                "\\end{figure}\n"
        clearpage = "\\clearpage\n"
        footer = "\\end{document}\n"
        self.pdf_dir = self.stat_dir+"pdf/"
        useful_funcs.create_dir(self.pdf_dir)
        f = open(self.pdf_dir+"main.tex","w")
        f.write(header)
        for i in [4,3,2]:
            quasars_times = quasars[quasars["yes"]==i]
            f.write(section % (str(i)+"bands, $>$99.7\%"))
            for index, row in quasars_times.iterrows():
                name = row["name"]
                useful_funcs.create_dir(self.pdf_dir+name)
                os.system("cp "+self.output_dir+name+"/periodogram.png "+\
                          self.pdf_dir+name)
                os.system("cp "+self.output_dir+name+"/lightcurve.png "+\
                          self.pdf_dir+name)
                os.system("cp "+self.output_dir+name+"/ACF.png "+\
                          self.pdf_dir+name)
                f.write(body1 % (name,name,name))
            f.write(clearpage)
        f.write(footer)
        f.close()

    ####################
    # Auto-correlation #
    ####################

    def _confidence_level(self,N,level=0.997):

        import scipy.special

        sigma = scipy.special.erfinv(1.-2.*(1-level))*np.sqrt(2)
        R_p = (-1+sigma*np.sqrt(N-2))/(N-1)
        R_n = (-1-sigma*np.sqrt(N-2))/(N-1)

        return R_p,R_n

    def _ACF_mock_confidence_level(self,power_mock):

        power_mock_at_each_period = zip(*power_mock)
        percentage_u = 99.7
        percentage_l = 68.0
        percentile_line_u = np.percentile(power_mock_at_each_period,\
                            percentage_u,axis=1)
        percentile_line_l = np.percentile(power_mock_at_each_period,\
                            percentage_l,axis=1)

        return percentile_line_u,percentile_line_l

    def run_ACF(self,name,band,using_mock=False):

        # read light curves and save to the right format
        name_band = "%s_%s" % (name,band)
        lc = np.load("%s/real/lc_%s.npy" % (self.output_dir+name,band))
        np.savetxt("%s.csv" % name_band, lc.T)

        # run ZDCF
        if using_mock:
            real_lc = np.load("%s/real/lc_%s.npy" % (self.output_dir+name,band))
            period, error = real_lc[0],real_lc[2]
            lcs = np.load("%s/mock/lcs_%s.npy" % (self.output_dir+name,band))
            power_all = []
            for lc in lcs: 
                mock_lc = np.array([period,lc,error]).T
                np.savetxt("%s.csv" % name_band,mock_lc)

                os.system("echo %s | run_ACF " % name_band)
                os.system("rm %s.lc1 %s.lc2 %s.csv" % \
                          (name_band,name_band,name_band) )
                data = np.genfromtxt("%s.dcf" % name_band,dtype=float)
                power = data[:,3]
                os.system("rm %s.dcf" % name_band)

                power_all.append(power)
            np.save("%s/mock/acf_%s.npy" % \
                    (self.output_dir+name,band),power_all)
            
        else:
            lc = np.load("%s/real/lc_%s.npy" % (self.output_dir+name,band))
            np.savetxt("%s.csv" % name_band,lc.T)

            os.system("echo %s | run_ACF " % name_band)
            os.system("rm %s.lc1 %s.lc2 %s.csv" % \
                     (name_band,name_band,name_band) )
            data = np.genfromtxt("%s.dcf" % name_band,dtype=float)
            os.system("rm %s.dcf" % name_band)
            np.save("%s/real/acf_%s.npy" % \
                    (self.output_dir+name,band),data)


    def generate_ACF_results(self,quasar):

        name = quasar["name"]

        # read lightcurves and save to the right format
        for band in self.band_list: 
            self.run_ACF(name,band)
            self.run_ACF(name,band,using_mock=True)

        # plot ACF result
        
        ACF = plot(4,1,figsize=(8,8),sharex=True)
        for band in self.band_list:
            if os.path.exists("%s/real/acf_%s.npy" % \
                              (self.output_dir+name,band)):
                data = np.load("%s/real/acf_%s.npy" % \
                               (self.output_dir+name,band))
                mock_data = np.load("%s/mock/acf_%s.npy" % \
                                    (self.output_dir+name,band))
                time = data[:,0]/365.
                xerr_l = data[:,1]/365.
                xerr_u = data[:,2]/365.
                power = data[:,3]
                yerr_l = data[:,4]
                yerr_u = data[:,5]
                N = data[:,6]
                #boundary_u,boundary_l =  self._confidence_level(N,0.997)
                boundary_u,boundary_l = self._ACF_mock_confidence_level(\
                                        mock_data)
                print boundary_u,boundary_l
                popt,perr,xn,yn,r_squared = self._fitting_exp_sin(time,power)

                ACF.plot_ACF(time,xerr_l,xerr_u,power,yerr_l,yerr_u,\
                         boundary_u,boundary_l,band)
                ACF.plot_fit_curve(xn,yn,band)
                ACF.plot_year_estimate(popt,perr,r_squared,band)

        ACF.savefig(self.output_dir+name,"/ACF.png",name)

        '''
        from plot import plot_ACF
        from matplotlib import pyplot as plt
        f,ax = plt.subplots(1,1,figsize=(12,3))
        for band in band_list:
            name_band = "%s_%s" % (name,band)
            if os.path.exists("%s.dcf" % name_band):
                data = np.genfromtxt("%s.dcf" % name_band,dtype=float)
                os.system("rm %s.dcf" % name_band)
                time = data[:,0]/365.
                xerr_l = data[:,1]/365.
                xerr_u = data[:,2]/365.
                power = data[:,3]
                yerr_l = data[:,4]
                yerr_u = data[:,5]
                N = data[:,6]
                boundary_u,boundary_l =  self._confidence_level(N,0.997)

                plot_ACF(ax,time,xerr_l,xerr_u,power,yerr_l,yerr_u,\
                         boundary_u,boundary_l,band)
        ax.set_xlabel("time(yr)")
        ax.set_ylabel("ACF")
        ax.legend()
        f.tight_layout(rect=[0, 0.03, 1, 0.95])
        f.suptitle(name)
        f.savefig(self.output_dir+name+"/ACF.png",dpi=300)
        '''

