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

        self.f,self.axes  = plt.subplots(cols, rows,**args)
        if cols==4 and rows ==1:
            self.ax_list = {"g":self.axes[0],"r":self.axes[1],\
                            "i":self.axes[2],"z":self.axes[3]}
        elif cols==2 and rows ==2:
            self.ax_list = {"g":self.axes[0,0],"r":self.axes[0,1],\
                            "i":self.axes[1,0],"z":self.axes[1,1]}
        self.color_list = {"g":"g","r":"orange",\
                           "i":"brown","z":"purple"}
        self.fmt_list = {"DES":{"fmt":"o","markersize":5},\
                         "SDSS_corr":{"fmt":"s","markersize":5},\
                         "PS":{"fmt":"^","markersize":5}
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
        ax.fill_betweenx([0.0, 1.05], max(_freq)/365/3,  max(_freq)/365, color='grey', alpha='0.5')
        if band == "i" or  band == "z":
            ax.set_xlabel("Period(yr)")
        if band == "g" or  band == "i":
            ax.set_ylabel("Power")
        ax.annotate(band, xy=(0.95, 0.90),xycoords='axes fraction',\
                    size=12, ha='right', va='top', color=color_list[band],
                    bbox=dict(boxstyle='round', fc='w'))
#        ax.legend()


    def plot_mock_periodogram(self,_freq, psd,band):

        ax_list = {"g":self.axes[0,0],"r":self.axes[0,1],\
                   "i":self.axes[1,0],"z":self.axes[1,1]}
        color_list = {"g":"g","r":"orange",\
                      "i":"brown","z":"purple"}
        ax = ax_list[band]
        ax.plot(_freq/365,psd,label=band,c="grey",linewidth=0.01)

    def plot_boost_periodogram(self,_freq, psd,error,band):

        sigma_level = 3
        upper = psd+error*sigma_level
        lower = psd-error*sigma_level
        ax = self.ax_list[band]
        ax.plot(_freq/365,upper,label=band,c=self.color_list[band],\
                linewidth=0.5)
        ax.plot(_freq/365,lower,label=band,c=self.color_list[band],\
                linewidth=0.5)

    def plot_peak_period(self,period):

        ax.axvline(x=period, color="r")


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
           ax.plot(_freq/365,bounday_psd_at_each__freq,"--",c="black",linewidth=0.5)

    def plot_light_curve(self,time,signal,error,survey,band):

        if survey in self.fmt_list.keys() : fmt = self.fmt_list[survey]
        else : fmt = {"fmt":"x","markersize":5}
        ax = self.ax_list[band]
        ax.errorbar(time,signal,yerr=error,label=band,c=self.color_list[band],\
                    **fmt)
        bottom,top = ax.get_ylim()
        if bottom > np.min(signal)-0.5: bottom = np.min(signal)-0.5
        if top < np.max(signal)+0.5: top = np.max(signal)+0.5
        ax.set_ylim(bottom,top)
        print survey,min(time),max(signal)
        if band == "z":
            ax.set_xlabel("MJD")
        #ax.set_ylabel("Mag") 
        ax.set_ylabel("Flux(nanomaggy)")
        ax.annotate(band, xy=(0.98, 0.9),xycoords='axes fraction',\
                    size=12, ha='right', va='top', color=self.color_list[band],\
                    bbox=dict(boxstyle='round', fc='w'))

    def plot_fit_curve(self,time,signal,band):

        ax = self.ax_list[band]
        ax.plot(time,signal,label=band,c=self.color_list[band],linestyle="--",\
                linewidth=1)

    def plot_walkers(self,samples):

        labels = ["tau", "c", "b"]
        for i in range(2):
            ax = self.axes[i]
            ax.plot(samples[:, :, i].T, "k", alpha=0.3)
            ax.set_ylabel(labels[i])
            #ax.yaxis.set_label_coords(-0.1, 0.5)
        self.axes[-1].set_xlabel("step number")
        

    def plot_mock_curve(self,time,signal,band):

        ax_list = {"g":self.axes[0],"r":self.axes[1],\
                   "i":self.axes[2],"z":self.axes[3]}
        ax = ax_list[band]
        ax.plot(time,signal,label=band,c="grey",linewidth=0.001)

    def savefig(self,dir_output,name,title):

        self.f.tight_layout(rect=[0, 0.03, 1, 0.95])
        self.f.suptitle(title)
        self.f.savefig(dir_output+name)
        plt.close(self.f)


        

