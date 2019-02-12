import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('agg')
from javelin.zylc import get_data
from javelin.lcmodel import Cont_Model
from quasar_drw import quasar_drw as qso_drw
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
        self.random_state = np.random.RandomState(0)
        self.band_list = ["g", "r", "i", "z"]
        self.surveys = ["DES", "SDSS_corr"]
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
        np.savetxt(filename,zipdata,delimiter=" ",comments="")

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
        data = pd.read_csv(file_path, comment="#")
        return data

    def error_boostraping(self, lc, periodogram, band):

        if self.test is True: number_times = 100
        else: number_times = 10000
        period_obs, signal_obs = lc.periodogram(lc.time, lc.signal)
        signal_boost = self.random_state.normal(lc.signal, lc.error, \
                       [number_times, len(lc.signal)])
        signal_total = []
        for signal in signal_boost:
            period_boost, psd_boost = lc.periodogram(lc.time, signal)
            signal_total.append(psd_boost)

        signal_total = np.array(signal_total, dtype=np.float64)

        error_list = []
        for i in range(len(signal_total[0])):
            error = np.std(signal_total[:, i])
            error_list.append(error)
        error_list = np.array(error_list, dtype=np.float64)
        periodogram.plot_boost_periodogram(period_obs, signal_obs, error, band)
        
        return error_list

    def check_period_max_amp(self, period, psd):
        # only search for peak amp between 500 days and 3 cycles.
        interval = (period<max(period)/3.) & (period>self.period_lowerlim)
        period_within = period[interval]
        amp_within = psd[interval]
        period_max = period_within[np.where(amp_within==np.max(amp_within))]
        return period_max


    def fitting(self, time, signal, error, period):

        def sin_func(x, amplitude, ref_day, median):
            return amplitude*np.sin(2*np.pi*x/period+ref_day)+median
        p0 = [1, 50000, 20]
        popt, pcov = curve_fit(sin_func, time, signal, p0=p0)

        xn = np.linspace(np.min(time)-100, np.max(time)+100, 10000)
        yn = sin_func(xn, *popt)
        return xn, yn

    def clean_parameters_list(self, parameters_list):

        num_parameters = len(parameters_list[0])
        excluded = np.zeros(len(parameters_list[:,0]), dtype=bool)
        for i in range(0, num_parameters-1):
            array = parameters_list[:, i]
            upper = np.percentile(array, 90)
            lower = np.percentile(array, 10)
            excluded = (excluded) | ( ~((array<upper) & (array>lower) ) ) 
        parameters_list = parameters_list[~excluded]
        return parameters_list

    def save_period_amp(self, _freq, psd, error, filename):

        save_data = zip(_freq, psd, error)
        np.savetxt(filename, save_data, delimiter=",", comments="",\
                   header="period,amplitude,error")

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

    def tailored_simulation(self,lc,band,z,name,periodogram,lightcurve):

        psd_mock_all = []
        if self.test is True: 
            nwalkers,burnin,Nsteps,draw_times = 100,20,100,10
        else:
            nwalkers = 500
            burnin = 150
            Nsteps = 500
            draw_times = 5000

        samples =  lc.fit_drw_emcee(nwalkers=nwalkers, burnin=burnin,\
                           Nstep=Nsteps,random_state=self.random_state)
        walkers = plot(2,1,figsize=(10,7), sharex=True)
        walkers.plot_walkers(np.exp(samples))
        walkers.savefig(self.output_dir+name,"/walker_"+band+".png",band+" band")

        parameters_list = samples[:, burnin:, :].reshape((-1, 2))
        parameters_list_good = self.clean_parameters_list(parameters_list)

        for i in range(draw_times):
            tau,c = np.exp(parameters_list_good[self.random_state.randint(\
                    len(parameters_list_good))])
            mock_time,mock_signal = lc.generate_mock_lightcurve(tau,c,lc.time,\
                                    lc.signal,z,random_state=self.random_state)
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

    def tailored_simulation_javelin(self,lc,band,z,name,periodogram,lightcurve):

        useful_funcs.create_dir(self.lc_dir+"combined")
        useful_funcs.create_dir(self.lc_dir+"combined/"+name)

        filename = self.make_combined_lc_for_javelin(lc,name,band)
        javdata = get_data(filename,names=["Continuum"])
        cont = Cont_Model(javdata)
        cont.do_mcmc(nwalkers=100, nburn=50, nchain=100,fchain="test.dat")


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
                data = self.add_flux_fluxerr_into_data(data)
                lc = qso_drw(data["mjd_obs"],data["flux_psf"],\
                             data["flux_err_psf"], quasar["z"],\
                             preprocess=True)
                time, signal, error = lc.get_lc()

                lightcurve.plot_light_curve(time,signal,error,survey,band)

                # add other surveys
                for survey in self.surveys[1:]:
                    data = self.read_lightcurve(name,survey,band)
                    data = self.add_flux_fluxerr_into_data(data)
                    lc2 = qso_drw(data["mjd_obs"],data["flux_psf"],\
                                 data["flux_err_psf"], quasar["z"],\
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
                    error_list = self.error_boostraping(lc,periodogram,band)

                    period_max = self.check_period_max_amp(period, psd)
                    if period_max is  None:
                        print ("Can not find the peak value !!")
                    else:
                        xn,yn = self.fitting(lc.time,lc.signal,lc.error,period_max)
                        lightcurve.plot_fit_curve(xn,yn,band)

                    self.tailored_simulation(lc,band,quasar["z"],name,\
                                             periodogram,lightcurve)
                    self.save_period_amp(period, psd,error_list,\
                                         self.output_dir+name+\
                                         "/periodogram_"+band+".csv")
        lightcurve.savefig(self.output_dir+name,"/lightcurve.png",name)
        periodogram.savefig(self.output_dir+name,"/periodogram.png",name)


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
                if row["peak_"+band]<1:
                    yes= yes+1
            cand.at[index,"yes"] = yes
        strong_cand = cand[cand["yes"]>1]
        strong_cand.to_csv(self.stat_dir+"strong_candidates.csv",index=False)

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
                "\\end{figure}\n"
        clearpage = "\\clearpage\n"
        footer = "\\end{document}\n"
        self.pdf_dir = self.stat_dir+"pdf/"
        useful_funcs.create_dir(self.pdf_dir)
        f = open(self.pdf_dir+"main.tex","w")
        f.write(header)
        for i in [4,3,2]:
            quasars_times = quasars[quasars["yes"]==i]
            f.write(section % (str(i)+"bands, $>$99\%"))
            for index, row in quasars_times.iterrows():
                name = row["name"]
                useful_funcs.create_dir(self.pdf_dir+name)
                os.system("cp "+self.output_dir+name+"/periodogram.png "+\
                          self.pdf_dir+name)
                os.system("cp "+self.output_dir+name+"/lightcurve.png "+\
                          self.pdf_dir+name)
                f.write(body1 % (name,name))
            f.write(clearpage)
        f.write(footer)
        f.close()
