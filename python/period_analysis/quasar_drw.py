import numpy as np
from numpy import mean, median
from astropy.stats import sigma_clip
import scipy.signal as sig
#from lomb_scargle_red_fix import lomb
from astroML.time_series import lomb_scargle
import os, subprocess
from scipy.optimize import curve_fit
import emcee
import scipy.optimize as op
from scipy.signal import medfilt 
from astroML.time_series import generate_damped_RW
from scipy.optimize import minimize, rosen, rosen_der
from plot import plot
#import multiprocessing as mp
import pyLCSIM


class quasar_drw:

    """ 
    Originally written by Wei-Ting Liao 2018
    Modified by Yu-Ching (Tony) Chen Oct. 2018
    This version directly fits state space 
    """
    
    def __init__(self, time, signal, error, redshift, preprocess=True):
        self.time     = np.array(time,   dtype=np.float64)
        self.signal   = np.array(signal, dtype=np.float64)
        self.error    = np.array(error,  dtype=np.float64)
        self.redshift = float(redshift)
        self._preprocessed = False

        if ( len(time) != len(signal) ) or ( len(time)!= len(error) ):
            print("[quasar_drw] Error in input data: time, signal, error must have the same length.")
        
        self._sort_data()
        if preprocess == True:
            self._preprocess()
        self.__initiate()

    
    def __initiate(self):
        ## parameters for periodogram
        if (len(self.time) >= 2 and float( np.max(self.time) - np.min(self.time) ) > 0.0 ):
            self.__Tspan    = float( np.max(self.time) - np.min(self.time) )
            self.__Ndata    = len(self.signal)
            self.__psd_freq = \
                              np.linspace(1.0/self.__Tspan, 1./(0.75*365), self.__Ndata)
               # np.linspace(1.0/self.__Tspan, self.__Ndata/(2.0*self.__Tspan), 10*self.__Ndata)
               # np.linspace(1.0/self.__Tspan, self.__Ndata/(2.0*self.__Tspan), self.__Ndata) 

            self.__dt = self.__Tspan / float(self.__Ndata)
            self.__df = self.__psd_freq[1] - self.__psd_freq[0]
        else:
            pass
        
    
    def get_Tspan(self):
        if (len(self.time) >= 2):
            return self.__Tspan
        else:
            pass
            
    def get_Ndata(self):
        if (len(self.time) >= 2):
            return self.__Ndata
        else:
            pass
            
    def get_psd_freq(self):
        return self.__psd_freq
        
    def get_psd_time(self):
        return 1.0/self.__psd_freq
        
    def get_psd_time_err(self):
        period = self.get_psd_time()
        return (self.__df * period**2.)/2.
            
    
    def _preprocess(self):
        #self._sort_data()
        self._no_outlier()
        self._bin_data()
        self._preprocessed = True
        self.__initiate()
        
    
    
    def get_lc(self):
        """ output: time, signal, error """
        return (self.time, self.signal, self.error)
        
        
    def get_redshift(self):
        return self.redshift
        
    
    def add_lc(self, time, signal, error, preprocess=True):
        self.time   = np.array(np.append(self.time,   time),   dtype=np.float64) 
        self.signal = np.array(np.append(self.signal, signal), dtype=np.float64)
        self.error  = np.array(np.append(self.error,  error),  dtype=np.float64)
        
        self._sort_data()
        
        self._preprocessed = False
        
        if (preprocess == True):
            self._preprocess()
            
        self.__initiate()
        
         
        
    def ls_astroML(self):
        """
        calculate periodogram using generalized Lomb-Scargle periodogram from AstroML
        function description: http://www.astroml.org/modules/generated/astroML.time_series.lomb_scargle.html
        example: http://www.astroml.org/book_figures/chapter10/fig_LS_example.html
        """
        LS_lc = lomb_scargle(self.time, self.signal, self.error, self.__psd_freq*(2.0*np.pi), generalized=True)
        
        return 1.0/self.__psd_freq, LS_lc

    def periodogram(self,time,signal):

        LS_lc = lomb_scargle(time, signal, self.error, self.__psd_freq*(2.0*np.pi), generalized=True)
        return 1.0/self.__psd_freq, LS_lc

    ### ********************************* ###
    ###         for drw modeling          ###
    ### ********************************* ###
    
    def fit_drw_emcee(self, nwalkers=500, burnin=150, Nstep=500,random_state=np.random.RandomState(0)):
        #ndim = 2
        ndim    = 3
        pos     = []
        
        z           = self.redshift
        time        = self.time
        signal      = self.signal
        error       = self.error
        
        # use most likely val as a initial guess
        nll = lambda *args: -lnlike(*args)
#        signal = signal-np.mean(signal)
        result = op.minimize(nll, [np.log(300.), np.log(0.001), np.log(np.mean(signal)/300.)], args=(self.time, self.signal, self.error, self.redshift),method="Nelder-Mead")
        # "BFGS"

        tau_center = np.exp(result["x"][0])
        c_center   = np.exp(result["x"][1])
        b_center   = np.exp(result["x"][2])
        #c_center = 0.0001

        
        print("Initial guess of (tau, c, b) = (" + format(np.exp(result["x"][0]), ".2f") + ", " \
                                                 + format(np.exp(result["x"][1]), ".2e") + ", " \
         #                                        )
                                                 + format(np.exp(result["x"][2]), ".2f") + " )" )
        
        ## initiate a gaussian distribution aroun dthe mean value
        ## modify this part if needed
#        tau_center = 300
#        c_center = 0.02
        tau_sample = random_state.lognormal(mean=np.log(tau_center), sigma=1, size=nwalkers)
#        tau_sample = np.random.lognormal(mean=np.log(tau_center), sigma=0.1, size=nwalkers)
        c_sample   = random_state.lognormal(mean=np.log(c_center),   sigma=1, size=nwalkers)
#        c_sample   = np.random.lognormal(mean=np.log(c_center),   sigma=0.1, size=nwalkers)
        b_sample   = random_state.normal(tau_center*b_center, 0.1,\
                     size=nwalkers)/tau_sample
        #random_state.lognormal(mean=np.log(b_center),   sigma=0.2, size=nwalkers)
#        b_sample   = np.random.lognormal(mean=np.log(b_center),   sigma=0.1, size=nwalkers)

        
        tau_sample, c_sample, b_sample = np.log(tau_sample), np.log(c_sample), np.log(b_sample)
#        tau_sample,c_sample = np.log(tau_sample), np.log(c_sample)
#        print list(tau_sample)
        
        for i in range(nwalkers):
            parameter = np.array([tau_sample[i], c_sample[i], b_sample[i]])
#            parameter = np.array([tau_sample[i], c_sample[i]])
            pos.append(parameter)
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(time, signal, error, z), a=2.0)
        
        # import random state 
        sampler.random_state = random_state.get_state()
        # start MCMC
        sampler.run_mcmc(pos, Nstep)
    
        # remove burn-in
        burnin = burnin
        #samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))

        ## depending on the preference, return whatever you prefer
        return sampler.chain

#        return [[tau_center,c_center,b_center]]

    def fit_drw_model_mcmc(self, nwalkers=500, burnin=150, Nstep=500,random_state=np.random.RandomState(0)):

        ndim    = 3
        pos     = []

        z           = self.redshift
        time        = self.time
        signal      = self.signal
        error       = self.error

        # use most likely val as a initial guess
        nll = lambda *args: -lnlike_drw(*args)
        result = op.minimize(nll, [np.log(300.), np.log(0.1), np.log(np.mean(signal))], args=(self.time, self.signal, self.error, self.redshift),method="Nelder-Mead")
        tau_center = np.exp(result["x"][0])
        sigma_center   = np.exp(result["x"][1])
        mean_center   = np.exp(result["x"][2])

        print("Initial guess of (tau, sigma, mean) = "+\
              "( %.2f, %.2f, %.2f )" % tuple(np.exp(result["x"]))
              )

       ## initiate a gaussian distribution aroun dthe mean value
        ## modify this part if needed
        pos = random_state.normal(loc=result["x"], scale=[0.1,0.1,0.01],\
                                  size=[nwalkers,len(result["x"])])
        #tau_sample = random_state.normal(loc=result["x"][0], scale=1, size=nwalkers)
        #sigma_sample   = random_state.normal(loc=result["x"][1],   scale=1, size=nwalkers)
        #mean_sample   = random_state.normal(loc=result["x"][2], scale= 0.1, size=nwalkers)

        #for i in range(nwalkers):
        #    parameter = np.array([tau_sample[i], sigma_sample[i], mean_sample[i]])
        #    pos.append(parameter)
        #pool = mp.Pool(mp.cpu_count()/6)
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_drw, args=(time, signal, error, z), a=2.0)#,pool=pool)

        # import random state 
        sampler.random_state = random_state.get_state()
        # start MCMC
        sampler.run_mcmc(pos, Nstep)
    
        # remove burn-in
        burnin = burnin

        #pool.close()
        #pool.join()

        ## depending on the preference, return whatever you prefer
        return sampler.get_chain(discard=burnin),sampler.get_log_prob(discard=burnin)

    def fit_periodic_model_mcmc(self, nwalkers=500, burnin=150, Nstep=500,random_state=np.random.RandomState(0),model="sin",name=None,band="g",drw_periodic=False):

        if drw_periodic:
            ndim    = 6
        else:
            ndim    = 5

        pos     = []

        z           = self.redshift
        time        = self.time - min(self.time)
        signal      = self.signal
        error       = self.error

        #q_data_path = "accretion"
        q_data_path = "."

        if model == "sin":
            self.sim_time = np.linspace(-1,10,100000)
            self.sim_signal = np.sin(self.sim_time*2*np.pi)
        elif model == "q011" or model == "q043":
            data1 = np.genfromtxt("%s/%s_primary.csv" % (q_data_path,model),dtype=float,delimiter=",")
            data2 = np.genfromtxt("%s/%s_secondary.csv" % (q_data_path,model),dtype=float,delimiter=",")
            self.sim_time = data1[:,0]
            if model == "q043": self.sim_time = self.sim_time/5.

            self.sim_time = self.sim_time-np.min(self.sim_time)-1
            self.sim_signal = data1[:,1]+data2[:,1]
            self.sim_signal = self.sim_signal - np.mean(self.sim_signal)
            

        #ini_pars = [np]

        # use most likely val as a initial guess
        nll = lambda *args: -lnprob_periodic(*args)
        if drw_periodic:
            result = op.minimize(nll, [np.log(1000),np.log(0.5),np.log(0.3),np.log(np.mean(signal)), np.log(300), np.log(0.3)], args=(time, signal, error, self.sim_time, self.sim_signal, self.redshift,drw_periodic),method="Nelder-Mead")
            #tau_center = np.exp(result["x"][0])
            #sigma_center   = np.exp(result["x"][1])
            #mean_center   = np.exp(result["x"][2])
            if result["x"][4] >  np.log(2000): result["x"][4] = np.log(300)
            if result["x"][5] < np.log(0.1): result["x"][5] = np.log(0.3)

            """
            if name == "J024613.89-004028.2": result["x"] = np.log(np.array([1170,52040,0.31,np.mean(signal),np.exp(2.55),np.exp(-1.7)]))
            if name == "J024703.24-010032.0": result["x"] = np.log(np.array([1620,54540,0.85,np.mean(signal),np.exp(2.85),np.exp(-0.4)]))
            if name == "J024944.66-000036.8": result["x"] = np.log(np.array([1100,54230,0.6,np.mean(signal),np.exp(2.58),np.exp(-0.9)]))
            if name == "J025406.26+002753.7": result["x"] = np.log(np.array([1490,53130,0.65,np.mean(signal),np.exp(2.9),np.exp(-0.5)]))
            """
            if name == "J024613.89-004028.2": result["x"] = np.log(np.array([1170,0.5,0.31,np.mean(signal),np.exp(2.55),np.exp(-1.7)]))
            if name == "J024703.24-010032.0": 
                result["x"] = np.log(np.array([1780,1.1,0.7,np.mean(signal),np.exp(2.5),np.exp(-0.7)]))
                if band == "r": result["x"][[2,4,5]] = [np.log(0.75),2.85,-0.5]
                elif band == "i": result["x"][[2,4,5]] = [np.log(0.6),3.1,-0.35]
                elif band == "z": result["x"][[2,4,5]] = [np.log(0.58),3.3,-0.3]
                if model == "q011": result["x"][[1,2]] = [np.log(1.42),np.log(0.8)]

            if name == "J024944.66-000036.8": 
                result["x"] = np.log(np.array([1150,1.25,0.6,np.mean(signal),np.exp(2.4),np.exp(-1.0)]))
                if model == "q011": result["x"][1] = np.log(0.9)
            if name == "J025406.26+002753.7": 
                result["x"] = np.log(np.array([1490,0.62,0.55,np.mean(signal),np.exp(2.82),np.exp(-0.61)]))
                if band == "r": result["x"][[2,4,5]] = [np.log(0.57),3.02,-0.45]
                elif band == "i": result["x"][[2,4,5]] = [np.log(0.55),3.05,-0.45]
                elif band == "z": result["x"][[2,4,5]] = [np.log(0.7),2.9,-0.65]
    
                if model == "q011": result["x"][1] = np.log(1.1)
            print("Initial guess of (t_ratio, t_shift, s_ratio, s_shift, "+\
                 "tau, sigma) = ( %.2f, %.2e, %.2f, %.2f, %.2f, %.2f )" % \
                 tuple(np.exp(result["x"])))

       ## initiate a gaussian distribution aroun dthe mean value
        ## modify this part if needed
            pos = random_state.normal(loc=result["x"], scale=[0.0001,0.0001,0.01,0.001,0.01,0.01],size=[nwalkers,len(result["x"])]) 
        else:
            result = op.minimize(nll, [np.log(1500),np.log(1.0),np.log(0.3),np.log(np.mean(signal)),np.log(0.05)], args=(time, signal, error, self.sim_time, self.sim_signal, self.redshift),method="Nelder-Mead")
            ## modify this part if needed
            if name == "J024613.89-004028.2":
                result["x"][[0,1]] = np.log(np.array([1170,0.8]))
                if model == "q011": result["x"][1] = np.log(1.14)
            if name == "J024703.24-010032.0":
                result["x"][[0,1]] = np.log(np.array([1780,0.9]))
                if model == "q011": result["x"][1] = np.log(1.25)
            if name == "J024944.66-000036.8": 
                result["x"][[0,1]] = np.log(np.array([1150,0.48]))
                if model == "q011": result["x"][1] = np.log(0.78)
            if name == "J025406.26+002753.7":
                result["x"][[0,1]] = np.log(np.array([1490,0.71]))
                if model == "q011": result["x"][1] = np.log(1.01)

            print("Initial guess of (t_ratio, t_shift, s_ratio, s_shift, add_error) = (%.2f, %.2e, %.2f, %.2f, %.2f )" % tuple(np.exp(result["x"])))
            pos = random_state.normal(loc=result["x"], scale=[0.01,0.01,0.01,0.001,0.001],size=[nwalkers,len(result["x"])])


        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_periodic, args=(time, signal, error, self.sim_time, self.sim_signal, z, drw_periodic), a=2.0)#, pool=pool)

        # import random state 
        sampler.random_state = random_state.get_state()
        # start MCMC
        sampler.run_mcmc(pos, Nstep)
    
        # remove burn-in
        burnin = burnin

        #pool.close()
        #pool.join()

        ## depending on the preference, return whatever you prefer
        return sampler.get_chain(discard=burnin),sampler.get_log_prob(discard=burnin)




    def generate_mock_lightcurve(self,tau,c,time,signal,z,random_state=np.random.RandomState(0)):

  
        time_res = time/(1+z)
        time_res_cont = np.linspace(min(time_res),max(time_res),int(max(time_res)-min(time_res)))
        xmean = np.mean(signal)
        SFinf = np.sqrt(c*tau/2.)
        lightcurve_DRW_res_cont = generate_damped_RW(time_res_cont,tau,z,xmean=xmean,SFinf=SFinf,random_state=random_state)
        lightcurve_DRW_res = np.interp(time_res, time_res_cont, lightcurve_DRW_res_cont)
        lightcurve_DRW_obs = lightcurve_DRW_res
        #lightcurve_DRW_obs = lightcurve_DRW_res_cont
        #time = time_res_cont*(1+z)
        return time,lightcurve_DRW_obs


    def get_BPL_lc(self,tau,time,signal):

        sim = pyLCSIM.Simulation()
        sim.addModel(myFunc, [0,-2,-3,1./tau,0.1])
        #myFunc(np.linspace(0.0001,10,10000),[0,-2,-3,0.003,0.1])

        rate_src    = 20
        rate_bkg    = 0
        t_exp       = 25*365
        dt          = 0.1
        frms        = 0.01
        nbins = t_exp/dt
        # Run the simulation
        sim.run(dt, nbins, rate_src, rms=frms)

        t,flux =  sim.getLightCurve()

        t = t+min(time)
        flux_BPL = np.interp(time, t, flux)

        return time,flux_BPL

    
    ### ********************************* ###
    ###  helper functions for preprocess  ###
    ### ********************************* ###
        
    def _sort_data(self):
        
        # take away points w/o data
        idx = self.error > 0.
        time   = self.time[idx]
        signal = self.signal[idx]
        error  = self.error[idx]
        
        idx = self.time > 0.
        time   = self.time[idx]
        signal = self.signal[idx]
        error  = self.error[idx]
        
        
        # sort
        idx = np.argsort(time)
        time   = time[idx]
        signal = signal[idx]
        error  = error[idx]
        
        # restore data
        self.time   = time
        self.signal = signal
        self.error  = error
    
    
    def _mask_sigma_clip_moving_avg(self,signal,sigma_level=5):

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

        
    def _no_outlier(self, sigma=3, iters=100):

        #idx = ((np.abs(self.signal) < 100.) & (self.signal > 0.)) # for magnitude
        idx = (self.signal > 0.) # for flux

        self.time   = self.time[idx]
        self.signal = self.signal[idx]
        self.error  = self.error[idx]
        
        #after_clip = sigma_clip(self.signal, sigma=sigma, iters=iters, cenfunc=median, copy=True)
        #idx = ~(after_clip.mask)

        idx = self._mask_sigma_clip_moving_avg(self.signal,sigma)

        self.time   = self.time[idx]
        self.signal = self.signal[idx]
        self.error  = self.error[idx]

        
    def _bin_data(self):
        time2   = []
        signal2 = []
        error2  = []
        count   = 0
    
        while(count < len(self.time)):
            idx = ( np.floor(self.time) == np.floor(self.time[count]) )
            signal_temp = self.signal[idx]
            error_temp  = self.error[idx]
            nn          = len(signal_temp)
        
            signal_temp, error_temp = self.__mag2flux(signal_temp, error_temp)
            signal_temp, error_temp = self.__weighted_mean(signal_temp, error_temp)
            signal_temp, error_temp = self.__flux2mag(signal_temp, error_temp)
        
            time2.append( np.floor(self.time[count]) ) 
            signal2.append( signal_temp )
            error2.append( error_temp )
        
            count += nn
        
        self.time   = np.asarray(time2)
        self.signal = np.asarray(signal2)
        self.error  = np.asarray(error2)
    
    
    ## bin input signal_temp to just one data point
    def _bin(self, signal_temp, error_temp):
        if (len(signal_temp) > 1):
            signal_temp, error_temp = self.__mag2flux(signal_temp, error_temp)
            signal_temp, error_temp = self.__weighted_mean(signal_temp, error_temp)
            signal_temp, error_temp = self.__flux2mag(signal_temp, error_temp)
        
        return signal_temp, error_temp
        

    def __mag2flux(self, signal, error):
        flux = 10.**(-1.*signal/2.5)
        return 10.**(-1.*signal/2.5), np.abs( -flux*error*np.log(10.)/2.5 )
    
    
    def __flux2mag(self, signal, error):
        return -2.5*np.log10(signal), np.abs( -2.5* error/signal/np.log(10.))
        
    
    def __weighted_mean(self, signal, error):
        signal_mean = np.sum(signal/error**2.) / np.sum(1./error**2.) 
        error_mean  = np.sqrt( np.sum(error**2.) ) / np.sqrt( np.float(len(signal)) )
        return signal_mean, error_mean
    
    ### *********************************** ###
    ###  END of helper func for preprocess  ###
    ### *********************************** ###
    
    

    
    
    ##### ------------------------------- #####
    ##### --- END of quasar_drw class --- #####
    ##### ------------------------------- #####


#######################################################
## ======== Function for MCMC model fitting ======== ##
#######################################################
def likelihood_a(time, tau_fit):
    Ndata = len(time)
    a_array = np.zeros(Ndata, dtype=np.float64)
    
    for i in range(1, Ndata):
        a_array[i] = np.exp( -(time[i]-time[i-1])/tau_fit )
        
    return a_array


def likelihood_omega(time, error, tau_fit, c_fit):
    Ndata   = len(time)
    a_array = likelihood_a(time, tau_fit)
    omega   = np.zeros(Ndata, dtype=np.float64)
    
    omega[0] = 0.5*tau_fit*c_fit
    
    for i in range(1, Ndata):
        omega[i] = omega[0]*(1.0-a_array[i]**2.0) + \
                    a_array[i]**2.0*omega[i-1]*( 1.0- omega[i-1]/(omega[i-1]+error[i-1]**2.0) )
    
    return omega
    

def likelihood_X(time, signal, error, tau_fit, c_fit, b_fit, X0=0.0):
    # this is X head (fitted X) in Kelly+09
    Ndata = len(time)
    a_array = likelihood_a(time, tau_fit)
    omega   = likelihood_omega(time, error, tau_fit, c_fit)
    signal_0mean = signal - b_fit*tau_fit
    
    X_array = np.zeros(Ndata, dtype=np.float64)
    X_array[0] = X0
    
    for i in range(1, Ndata):
        X_array[i] = a_array[i]*X_array[i-1] + \
                     a_array[i]*omega[i-1]/(omega[i-1]+error[i-1]**2.0)*(signal_0mean[i-1] - X_array[i-1])
    
    return X_array


def likelihood_P(time, signal, error, tau_fit, c_fit, b_fit):
    Ndata = len(time)
    a_array = likelihood_a(time, tau_fit)
    omega   = likelihood_omega(time, error, tau_fit, c_fit)
    X_array = likelihood_X(time, signal, error, tau_fit, c_fit, b_fit)
    signal_0mean = signal - b_fit*tau_fit
    
    P_array = np.zeros(Ndata, dtype=np.float64)
    for i in range(Ndata):
        P_array[i] = 1.0/np.sqrt((2.0*np.pi)*(omega[i]+error[i]**2.0)) * \
                     np.exp(-0.5 * ( (X_array[i]-signal_0mean[i])**2.0/(omega[i]+error[i]**2.0) ) )
                     
    return P_array


# set up for likelihood function
#def lnlike(theta, time, signal, error, z):
def lnlike(theta, time,signal,error,z):
    #lntau, lnc = theta
    lntau, lnc, lnb = theta
    #b_fit = np.mean(signal)/np.exp(lntau)
    #tau_fit = np.exp(lntau) * (1.0+z)
    #c_fit   = np.exp(lnc) * np.sqrt(1.0+z)
    #c_fit   = np.exp(lnc) / (1.0+z) 
    #b_fit   = np.exp(lnb) / (1.0+z)
    tau_fit = np.exp(lntau)
    c_fit = np.exp(lnc)
    b_fit = np.exp(lnb)
    P_array = likelihood_P(time, signal, error, tau_fit, c_fit, b_fit)
    Prob    = np.prod(P_array)
        
    if np.isfinite(Prob) and Prob > 0.:
        return np.log(Prob)
    else:
        return -np.inf


# set up for prior 
#def lnprior(theta, z, time):    
def lnprior(theta, z, time):
    # prior is determined in the rest frame, no need to multiply (1+z)
    lntau, lnc, lnb = theta
    tau_fit, c_fit, b_fit = np.exp(lntau), np.exp(lnc), np.exp(lnb)
#    lntau, lnc = theta
#    tau_fit,c_fit = np.exp(lntau),np.exp(lnc)
    
#    if 1.0 < tau_fit*(1.0+z) < (np.max(time)-np.min(time)) and c_fit > 0.0 :
    if c_fit > 0.0 and 10 < b_fit*tau_fit < 30 and 1.0 < tau_fit < (np.max(time)-np.min(time)): # mag
#    if c_fit> 0.0 and b_fit*tau_fit > 0.0001 and 1.0 < tau_fit < (np.max(time)-np.min(time)):
        return 0.0
    else:
        return -np.inf


# set up posterior 
def lnprob(theta, time, signal, error, z):    
    lp = lnprior(theta, z, time)
    lk = lnlike(theta, time, signal, error, z)
    lnprob_out = lp + lnlike(theta, time, signal, error, z)
    
    if ( np.isfinite(lp) and np.isfinite(lk) ):
        return lp+lk
    else:
        return -np.inf

####################
# model comparison #
####################

# set up posterior
def lnprob_drw (theta, time, signal, error, z):
    lp = lnprior_drw(theta, z, time)
    lk = lnlike_drw(theta, time, signal, error, z)
    #lnprob_out = lp + lnlike(theta, time, signal, error, z)

    if ( np.isfinite(lp) and np.isfinite(lk) ):
        return lp+lk
    else:
        return -np.inf

def lnprior_drw(theta, z, time):
    # prior is determined in the rest frame, no need to multiply (1+z)
    lntau, lnsigma, lnmean = theta
    tau_fit, sigma_fit, mean_fit = np.exp(lntau), np.exp(lnsigma), np.exp(lnmean)
    if tau_fit > 0.0 and 0 < mean_fit < 10**3 and 1.0 < tau_fit < (np.max(time)-np.min(time))/3.: # mag
        return 0.0
    else:
        return -np.inf

def lnprob_periodic (theta, time, signal, error, sim_time, sim_signal, z, periodic_drw=False):
    if periodic_drw:
        lp = lnprior_periodic(theta, z, time)
        lk = lnlike_periodic(theta, time, signal, error, sim_time, sim_signal, z)
    else:
        lp = lnprior_pure_periodic(theta,z ,time )
        lk = lnlike_pure_periodic(theta, time, signal,error, sim_time,sim_signal,z)
    #lnprob_out = lp + lnlike(theta, time, signal, error, z)

    if ( np.isfinite(lp) and np.isfinite(lk) ):
        return lp+lk
    else:
        return -np.inf

def lnprior_periodic(theta, z, time):
    # prior is determined in the rest frame, no need to multiply (1+z)
    lnt_ratio, lnt_shift, lns_ratio, lns_shift, lndrw_tau, lndrw_sigma = theta
    t_ratio, t_shift, s_ratio, s_shift, drw_sigma, drw_tau = np.exp(lnt_ratio),\
        np.exp(lnt_shift), np.exp(lns_ratio), \
        np.exp(lns_shift), np.exp(lndrw_sigma), np.exp(lndrw_tau)

    if drw_tau > 0.0  and 1.0 < drw_tau < (np.max(time)-np.min(time))/3. and \
       0 < s_shift < 10**3 and 0.05 < s_ratio < 10**3 and 0.2 < t_shift < 2 and\
       500. < t_ratio < (np.max(time)-np.min(time))/3. : # mag
        return 0.0
    else:
        return -np.inf

def lnprior_pure_periodic(theta, z, time):
    # prior is determined in the rest frame, no need to multiply (1+z)
    lnt_ratio, lnt_shift, lns_ratio, lns_shift, lnadd_error = theta
    t_ratio, t_shift, s_ratio, s_shift, add_error  = np.exp(lnt_ratio),\
        np.exp(lnt_shift), np.exp(lns_ratio), np.exp(lns_shift), np.exp(lnadd_error)

    if 0 < s_shift < 10**3 and 0.05 < s_ratio < 10**3 and 0.2 < t_shift < 2 and\
       500. < t_ratio < (np.max(time)-np.min(time))/3. and 0 < add_error < 1 : # mag
        return 0.0
    else:
        return -np.inf


def lnlike_periodic(theta, fit_time, fit_signal, fit_error, sim_time, sim_signal, z):

    t_ratio, t_shift, s_ratio, s_shift, drw_tau, drw_sigma  = np.exp(theta)

    fit = np.interp( fit_time, (sim_time-t_shift)*t_ratio, (sim_signal*s_ratio)+s_shift )
    model = fit_signal - fit

    """
    cov_D     = np.zeros([len(fit_signal),  len(fit_signal)])
    cov_sigma = np.zeros([len(fit_signal),  len(fit_signal)])


    for i in range(len(fit_time)):
        for j in range(len(fit_time)):
            delta_t = np.abs(fit_time[i] - fit_time[j])
            cov_D[i,j] = drw_sigma**2.0 * np.exp(- delta_t / drw_tau)

            if (i==j):
                cov_sigma[i,j] = fit_error[i]**2.0
    cov        = cov_D + cov_sigma

    """
    fit_time_x,fit_time_y = np.meshgrid(fit_time,fit_time)
    delta_t = np.abs(fit_time_x-fit_time_y)
    cov = drw_sigma**2.0 * np.exp(- delta_t / drw_tau)
    cov[np.arange(len(fit_time)),np.arange(len(fit_time))] += fit_error**2.0


    #
    cov_inverse = np.linalg.inv(cov)
    chi_square  = np.dot( model, np.dot(cov_inverse, model) )

    (sign, logdet) = np.linalg.slogdet( cov )

    lnlikeli = -0.5*logdet - 0.5*chi_square


    return lnlikeli


def lnlike_pure_periodic(theta, fit_time, fit_signal, fit_error, sim_time, sim_signal, z):

    t_ratio, t_shift, s_ratio, s_shift, add_error  = np.exp(theta)

    fit = np.interp( fit_time, (sim_time-t_shift)*t_ratio, (sim_signal*s_ratio)+s_shift )
    model = fit_signal - fit

    #error_2 = fit_error**2 + add_error**2
    #lnlikeli =  -0.5 * np.sum(model ** 2 / error_2 + np.log(error_2))
    #return lnlikeli


    cov = np.zeros((len(fit_time),len(fit_time)))
    #cov += error**2.0
    cov[np.arange(len(fit_time)),np.arange(len(fit_time))] += fit_error**2.0 + add_error**2.0


    #
    cov_inverse = np.linalg.inv(cov)
    chi_square  = np.dot( model, np.dot(cov_inverse, model) )

    (sign, logdet) = np.linalg.slogdet( cov )

    lnlikeli = -0.5*logdet - 0.5*chi_square


    return lnlikeli


def lnlike_drw(theta, fit_time, fit_signal, fit_error, z):


    tau, sigma, mean_mag = np.exp(theta)

    model = fit_signal - mean_mag


    """
    cov_D     = np.zeros([len(fit_signal),  len(fit_signal)])
    cov_sigma = np.zeros([len(fit_signal),  len(fit_signal)])


    for i in range(len(fit_time)):
        for j in range(len(fit_time)):
            if (i!=j):
                delta_t = np.abs(fit_time[i] - fit_time[j])
                cov_D[i,j] = sigma**2.0 * np.exp(- delta_t / tau)

            if (i==j):
                cov_sigma[i,j] = fit_error[i]**2.0
    cov        = cov_D + cov_sigma
    """
    fit_time_x,fit_time_y = np.meshgrid(fit_time,fit_time)
    delta_t = np.abs(fit_time_x-fit_time_y)
    cov = sigma**2.0 * np.exp(- delta_t / tau)
    cov[np.arange(len(fit_time)),np.arange(len(fit_time))] += fit_error**2.0


    #
    cov_inverse = np.linalg.inv(cov)
    chi_square  = np.dot( model, np.dot(cov_inverse, model) )

    (sign, logdet) = np.linalg.slogdet( cov )

    lnlikeli = -0.5*logdet - 0.5*chi_square


    return lnlikeli


            
########################################
## ======== END MCMCfunction ======== ##
########################################

def myFunc(f, p):

    # p: parameters [f1,f2,slope1,slope2,slope3]

    A = 1.
    f1 = f[f<p[3]]
    f2 = f[ (f<p[4]) & (f>p[3]) ]
    f3 = f[f>p[4]]

    P1 = A*f1**p[0]
    P2 = P1[-1]*f2**p[1]/(f2[0]**p[1])
    P3 = P2[-1]*f3**p[2]/(f3[0]**p[2])

    P = np.concatenate([P1,P2,P3])

    return P

