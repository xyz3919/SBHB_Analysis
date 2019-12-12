import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import numpy as np

class plot:

    """
    All the useful plot functions.
    
    by Yu-Ching (Tony) Chen
    ycchen@illinois.edu
    """

    def __init__(self, cols, rows,**args):

        plt.rc('font', family='serif')
        self.f,self.axes  = plt.subplots(cols, rows,**args)
        if cols==4 and rows ==1:
            self.ax_list = {"g":self.axes[0],"r":self.axes[1],\
                            "i":self.axes[2],"z":self.axes[3]}
        elif cols==2 and rows ==2:
            self.ax_list = {"g":self.axes[0,0],"r":self.axes[0,1],\
                            "i":self.axes[1,0],"z":self.axes[1,1]}
        self.color_list = {"g":"g","r":"orange",\
                           "i":"brown","z":"purple","total":"black"}
        self.fmt_list = {"DES":{"fmt":"o","markersize":5},\
                         "SDSS_corr":{"fmt":"s","markersize":5},\
                         "PS":{"fmt":"^","markersize":5},\
                         "ZTF":{"fmt":"+","markersize":5}
                        }

    def plot(self,x,y,band,log=False):

        color_list = {"g":"g","r":"orange",\
                      "i":"brown","z":"purple"}
        self.axes.plot(x,y,label=band,c=color_list[band])
        if log: 
            self.axes.set_xscale("log")
            self.axes.set_yscale("log")
        self.axes.legend()



    def plot_periodogram(self,_freq, psd,band):

        ax_list = {"g":self.axes[0,0],"r":self.axes[0,1],\
                   "i":self.axes[1,0],"z":self.axes[1,1]}
        color_list = {"g":"g","r":"orange",\
                      "i":"brown","z":"purple"}
        ax = ax_list[band]
        ax.plot(_freq/365,psd,label=band,c=color_list[band])
        ax.set_xlim(0.8,10)
        ax.set_ylim(0,1)
        ax.set_xscale("log")
        ax.set_xticks([1,2,4,8,10])
        ax.set_xticklabels([1,2,4,8,10])
        ax.fill_betweenx([0.0, 1.05], 0.8,  500./365., color='grey', alpha='0.5')
        ax.fill_betweenx([0.0, 1.05], max(_freq)/365/3,  max(_freq)/365, \
                         color='grey', alpha='0.5')
        if band == "i" or  band == "z":
            ax.set_xlabel("Period(yr)")
        if band == "g" or  band == "i":
            ax.set_ylabel("Power")
        ax.annotate(band, xy=(0.95, 0.90),xycoords='axes fraction',\
                    size=12, ha='right', va='top', color=color_list[band],
                    bbox=dict(boxstyle='round', fc='w'))
#        ax.legend()

    def plot_multi_periodogram(self,_freq,psd,band):

        ax = self.axes
        ax.plot(_freq/365,psd,label=band,c=self.color_list[band])
        ax.set_xlim(0.8,10)
        ax.set_ylim(0,1)
        ax.set_xscale("log")
        ax.set_xticks([1,2,4,8,10])
        ax.set_xticklabels([1,2,4,8,10])
        ax.fill_betweenx([0.0, 1.05], 0.8,  500./365., color='grey', alpha='0.5')
        ax.fill_betweenx([0.0, 1.05], max(_freq)/365/3,  max(_freq)/365, \
                         color='grey', alpha='0.5')
        ax.set_xlabel("Period(yr)")
        ax.set_ylabel("Power")
        #ax.annotate(band, xy=(0.95, 0.90),xycoords='axes fraction',\
        #            size=12, ha='right', va='top', color=color_list[band],
        #            bbox=dict(boxstyle='round', fc='w'))
        ax.legend()

    def plot_mock_periodogram(self,_freq, psds,band):

        ax_list = {"g":self.axes[0,0],"r":self.axes[0,1],\
                   "i":self.axes[1,0],"z":self.axes[1,1]}
        color_list = {"g":"g","r":"orange",\
                      "i":"brown","z":"purple"}
        ax = ax_list[band]
        for psd in psds:
            ax.plot(_freq/365,psd,label=band,c="lightgrey",linewidth=0.001)

    def plot_boost_periodogram(self,_freq, psd,error,band):

        sigma_level = 3
        upper = psd+error*sigma_level
        lower = psd-error*sigma_level
        ax = self.ax_list[band]
        ax.plot(_freq/365,upper,label=band,c=self.color_list[band],\
                linewidth=0.5)
        ax.plot(_freq/365,lower,label=band,c=self.color_list[band],\
                linewidth=0.5)

    def plot_peak_period(self,period, sig_level, band):

        cutout = 0.26
        ax = self.ax_list[band]

        bool_array = np.array(sig_level)< cutout 
        period_within = (period > 500 ) & (period < max(period)/3.)
        bool_array = period_within & bool_array
        boundaries = period[np.append(bool_array[0],np.diff(bool_array))]

        if len(boundaries) % 2 == 1: boundaries=np.append(boundaries,period[-1])
        for i in range(len(boundaries)/2):
            ax.axvspan(boundaries[i*2+1]/365,boundaries[i*2]/365,color='blue',\
                       alpha = 0.5)

    def plot_confidence_level(self,_freq, psd_total,band):
        ax_list = {"g":self.axes[0,0],"r":self.axes[0,1],\
                   "i":self.axes[1,0],"z":self.axes[1,1]}
        color_list = {"g":"g","r":"orange",\
                      "i":"brown","z":"purple"}
        ax = ax_list[band]
        psd_at_each__freq = zip(*psd_total)         
        percentiles = [68.27,95.45,99.0,99.74,99.99]
        for percentile in percentiles:
           bounday_psd_at_each__freq = [np.percentile(psd,percentile) for psd in psd_at_each__freq]
           ax.plot(_freq/365,bounday_psd_at_each__freq,"--",c="black",linewidth=1)



    def plot_light_curve(self,time,signal,error,survey,band,adjust_lim=True,yaxis="mag"):

        if survey in self.fmt_list.keys() : fmt = self.fmt_list[survey]
        else : fmt = {"fmt":"x","markersize":5}
        ax = self.ax_list[band]
        ax.errorbar(time,signal,yerr=error,label=band,c=self.color_list[band],\
                    **fmt)
        #bottom,top = ax.get_ylim()
        #if bottom > np.min(signal)-0.1: bottom = np.min(signal)-0.1
        #if top < np.max(signal)+0.1: top = np.max(signal)+0.1
        #ax.set_ylim(bottom,top)
        if band == "z":
            ax.set_xlabel("MJD")
        if yaxis== "mag":
            ax.set_ylabel("Magnitude") 
            bottom,top = ax.get_ylim()
            if top > bottom: 
                ax.set_ylim(ax.get_ylim()[::-1])
                bottom,top = ax.get_ylim()
            if bottom < np.max(signal)+0.1: bottom = np.max(signal+0.1)
            if top > np.min(signal)-0.1: top = np.min(signal)-0.1
            if adjust_lim:
                ax.set_ylim(bottom,top)
        else:
            ax.set_ylabel("Flux(nanomaggy)")
            bottom,top = ax.get_ylim()
            if bottom > np.min(signal)-0.1: bottom = np.min(signal)-0.1
            if top < np.max(signal)+0.1: top = np.max(signal)+0.1
            ax.set_ylim(bottom,top)


        ax.annotate(band, xy=(0.98, 0.9),xycoords='axes fraction',\
                    size=12, ha='right', va='top', color=self.color_list[band],\
                    bbox=dict(boxstyle='round', fc='w'))

    def plot_fit_curve(self,time,signal,band):

        ax = self.ax_list[band]
        ax.plot(time,signal,label=band,c=self.color_list[band],linestyle="--",\
                linewidth=1)

    def plot_walkers(self,samples):

        samples[:,:,1] = samples[:,:,0]*samples[:,:,1]/2
        samples[:,:,2] = samples[:,:,0]*samples[:,:,2]
        labels = ["tau", "var", "mean"]
        for i in range(3):
            ax = self.axes[i]
            upper,lower = np.percentile(samples[:,:,i],99),np.percentile(samples[:,:,i],1)
            ax.plot(samples[:, :, i].T, "k", alpha=0.05)
            ax.set_ylim(lower,upper)
            ax.set_ylabel(labels[i])
            ax.set_yscale('log')
            #ax.yaxis.set_label_coords(-0.1, 0.5)
        self.axes[-1].set_xlabel("step number")


    def plot_mock_curve(self,time,signal,band):

        ax_list = {"g":self.axes[0],"r":self.axes[1],\
                   "i":self.axes[2],"z":self.axes[3]}
        ax = ax_list[band]
        ax.plot(time,signal,label=band,c="grey",linewidth=0.1)

    def plot_ACF(self,x,xerr_l,xerr_u,y,yerr_l,yerr_u,boundary_u,\
                 boundary_l,band):
        ax = self.ax_list[band]
        ax.errorbar(x,y,yerr=[yerr_l,yerr_u],xerr=[xerr_l,xerr_u],\
                fmt=".",capsize=0.4,color=self.color_list[band],\
                label=band,elinewidth=0.2)
        ax.plot(x,boundary_u,linestyle='dashed',color="grey")
        ax.plot(x,boundary_l,linestyle='dashed',color="grey")
        if band == "z":
            ax.set_xlabel("Offset(yr)")
        ax.set_ylabel("Power")
        ax.set_xscale("log")
        ax.set_ylim(-1,1)
        ax.set_xlim(0.8,10)
        ax.set_xticks([1,2,4,8,10])
        ax.set_xticklabels([1,2,4,8,10])

    def plot_year_estimate(self,popt,perr,rsquared,band):
        
        period = popt[0]
        period_err  = perr[0]
        ax = self.ax_list[band]
        word = r"T = %.1f +- %.2f (yr) R$^2$=%.3f" % \
               (period,period_err,rsquared)
        ax.annotate(word, xy=(0.02, 0.9),xycoords='axes fraction',\
                    size=10, ha='left', va='top', color=self.color_list[band],\
                    bbox=dict(boxstyle='round', fc='w'))


    def savefig(self,dir_output,name,title):

        self.f.tight_layout(pad=1.05,rect=[0, 0.03, 1, 0.95])
        self.f.suptitle(title)
        self.f.savefig(dir_output+name,dpi=200)
        plt.close(self.f)


def plot_posterior(samples,likelihood,band,save_path,combine=True,model_comp=False):

    if model_comp:
        samples[:,0:2] = np.log10(samples[:,0:2])
    elif combine:
        samples[:,1] = samples[:,1]*samples[:,0]/2
        samples[:,2] = samples[:,2]*samples[:,0]
    else:
        samples[:,1] = samples[:,1]**2
        samples[:,0:2] = np.log10(samples[:,0:2])
    #print list(samples[:,2])
    import corner
    ndim = 3
    axrange = []
    sigma = np.std(samples,axis=0)
    median = np.median(samples,axis=0)
    maxlike = samples[likelihood == np.max(likelihood)][0]
#    for i in range(ndim):
#        axrange.append((median[i]-3*sigma[i],median[i]+3*sigma[i]))
    fig = corner.corner(samples,labels=[r"$log(\tau[days])$",\
                        r"$log(var[mag])$",r"$mean[mag]$"],\
                        quantiles=[0.16, 0.5, 0.84],\
                        show_titles=True, title_kwargs={"fontsize": 12},\
                        plot_datapoints=False,levels=(1-np.exp(-0.5),))
                        #range = axrange)
    axes = np.array(fig.axes).reshape((ndim, ndim))
#    value = np.median(samples, axis=0)
    # for the diagonal histograms
    for i in range(ndim):
        ax = axes[i, i]
        ax.axvline(maxlike[i], color="g")
    # for the 2-D posterior contours
    for yi in range(ndim):
        for xi in range(yi):
            ax = axes[yi, xi]
            ax.axvline(maxlike[xi], color="g")
            ax.axhline(maxlike[yi], color="g")
            ax.plot(maxlike[xi], maxlike[yi], "sg")
    fig.suptitle(band)
    fig.savefig(save_path)

def plot_posterior_drw_periodic(samples,likelihood,band,save_path):

    samples[:,5] = samples[:,5]**2
    samples[:,4:6] = np.log10(samples[:,4:6])
    #print list(samples[:,2])
    import corner
    ndim = 6
    axrange = []
    sigma = np.std(samples,axis=0)
    median = np.median(samples,axis=0)
    maxlike = samples[likelihood == np.max(likelihood)][0]
#    for i in range(ndim):
#        axrange.append((median[i]-3*sigma[i],median[i]+3*sigma[i]))
    fig = corner.corner(samples,labels=[r"t_ratio", r"t_shift", \
                        r"s_ratio", r"s_shift",r"$log(\tau[days])$",\
                        r"$log(var[mag])$"],\
                        quantiles=[0.16, 0.5, 0.84],\
                        show_titles=True, title_kwargs={"fontsize": 12},\
                        plot_datapoints=False,levels=(1-np.exp(-0.5),))
                        #range = axrange)
    axes = np.array(fig.axes).reshape((ndim, ndim))
#    value = np.median(samples, axis=0)
    # for the diagonal histograms
    for i in range(ndim):
        ax = axes[i, i]
        ax.axvline(maxlike[i], color="g")
    # for the 2-D posterior contours
    for yi in range(ndim):
        for xi in range(yi):
            ax = axes[yi, xi]
            ax.axvline(maxlike[xi], color="g")
            ax.axhline(maxlike[yi], color="g")
            ax.plot(maxlike[xi], maxlike[yi], "sg")
    fig.suptitle(band)
    fig.savefig(save_path)


def plot_posterior_carma(samples,band,save_path):

    samples[:,0:2] = np.log10(samples[:,0:2])
    #print list(samples[:,2])
    import corner
    ndim = 3
    axrange = []
    sigma = np.std(samples,axis=0)
    median = np.median(samples,axis=0)
#    for i in range(ndim):
#        axrange.append((median[i]-3*sigma[i],median[i]+3*sigma[i]))
    fig = corner.corner(samples,labels=[r"$log(\tau[days])$",\
                        r"$log(var[mag])$",r"$mean[mag]$"],\
                        quantiles=[0.16, 0.5, 0.84],\
                        show_titles=True, title_kwargs={"fontsize": 12},\
                        plot_datapoints=False)
                        #range = axrange)
    axes = np.array(fig.axes).reshape((ndim, ndim))
#    value = np.median(samples, axis=0)
    # for the diagonal histograms
    for i in range(ndim):
        ax = axes[i, i]
        ax.axvline(median[i], color="g")
    # for the 2-D posterior contours
    for yi in range(ndim):
        for xi in range(yi):
            ax = axes[yi, xi]
            ax.axvline(median[xi], color="g")
            ax.axhline(median[yi], color="g")
            ax.plot(median[xi], median[yi], "sg")

    fig.suptitle(band)
    fig.savefig(save_path)



def plot_drw_parameters(tau,sigma,band,save_path):

    import corner
    #samples = np.array([sigma,tau]).T
    #fig = corner.corner(samples,labels=[r"$log(var[mag])$",\
    #      r"$log(\tau[day])$"],titles_kwargs={"fontsize":12})
    SF = sigma+np.log10(np.sqrt(2))
    samples = np.array([SF,tau]).T
    fig = corner.corner(samples,labels=[r"$log(SF(mag))$",\
          r"$log(\tau[days])$"],titles_kwargs={"fontsize":12},\
          range=[(-2.0,0.),(0.,4.)],\
          quantiles=[0.16, 0.5, 0.84],show_titles=True)
    fig.suptitle(band)
    fig.savefig(save_path)

def plot_ACF(ax,x,xerr_l,xerr_u,y,yerr_l,yerr_u,boundary_u,boundary_l,band):

    color_list = {"g":"g","r":"orange",\
                  "i":"brown","z":"purple","total":"black"}
    ax.errorbar(x,y,yerr=[yerr_l,yerr_u],xerr=[xerr_l,xerr_u],\
                fmt=".",capsize=0.5,color=color_list[band],\
                label=band,elinewidth=0.2)
    ax.plot(x,boundary_u,linestyle='dashed',color="grey")#color_list[band])
    ax.plot(x,boundary_l,linestyle='dashed',color="grey")#color_list[band])
    ax.set_ylim(-1,1)

        

