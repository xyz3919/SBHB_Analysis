import matplotlib
matplotlib.use('Agg')
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
        self.color_list = {"g":"g","r":"orange",\
                           "i":"brown","z":"purple"}
        if cols==4 and rows ==1:
            self.ax_list = {"g":self.axes[0],"r":self.axes[1],\
                            "i":self.axes[2],"z":self.axes[3]}
        if cols == 2 and rows == 2:
            self.ax_list = {"g":self.axes[0,0],"r":self.axes[0,1],\
                            "i":self.axes[1,0],"z":self.axes[1,1]}
        self.fmt_list = {"DES":{"fmt":"o","markersize":5},\
                         "SDSS_corr":{"fmt":"s","markersize":5},\
                         "PS":{"fmt":"^","markersize":5},\
                         "ZTF":{"fmt":"x","markersize":5}
                        }

    def plot(self,x,y,band,log=False):

        self.axes.plot(x,y,label=band,c=self.color_list[band])
        if log:
            self.axes.set_xscale("log")
            self.axes.set_yscale("log")
        self.axes.legend()

    def plot_histogram(self,number_array,band):

        ax = self.ax_list[band]
        ax.set_yscale("log",nonposy='clip')
        ax.hist(number_array,bins=40,\
                range=(0,np.percentile(number_array,99)*1.05),\
                facecolor = self.color_list[band])
        if band == "i" or  band == "z":
            ax.set_xlabel("Number of epochs",fontsize=14)
        if band == "g" or  band == "i":
            ax.set_ylabel("Number of quasars",fontsize=14)
        ax.set_xlim(0,np.percentile(number_array,99)*1.05)
        ax.text(0.95, 0.95, band,size=16,color=self.color_list[band],\
                horizontalalignment='right',verticalalignment='top',\
                transform = ax.transAxes)

               
        #ax.annotate(band, xy=(-12, -12), xycoords='axes points',\
        #            size=12, ha='right',va='top',color=self.color_list[band],\
        #            bbox=dict(boxstyle='round', fc='w'))
    def plot_filter(self,x,y,band,camera):

        if camera == "SDSS": linestyle = "--"
        elif camera == "DES": linestyle = "-" 
        self.axes.plot(x,y,label=camera+" "+band,linestyle=linestyle,\
                       c=self.color_list[band])
        self.axes.set_xlabel("Wavelength(A)")
        self.axes.set_ylabel("Truoghtout")
        self.axes.legend(prop={'size': 10})

    def plot_redshift_mag(self,redshift_raw,mag_raw):

        from mpl_toolkits.axes_grid1 import make_axes_locatable
        import matplotlib.ticker as plticker
        mag = mag_raw[((mag_raw < 50) & (mag_raw > -50))]
        redshift = redshift_raw[((mag_raw < 50) & (mag_raw > -50))]

        self.axes.scatter(redshift, mag,c="black",s=6)
        self.axes.set_xlabel("Redshift",fontsize=14)
        self.axes.set_ylabel("i Band Magnitude",fontsize=14)
        divider = make_axes_locatable(self.axes)
        self.axHistx = divider.append_axes("top", 1.4, pad=0.1, \
                       sharex=self.axes)
        self.axHisty = divider.append_axes("right", 1.4, pad=0.1, \
                       sharey=self.axes)
        # make some labels invisible
        plt.setp(self.axHistx.get_xticklabels() + \
                 self.axHisty.get_yticklabels(),visible=False)
        self.axes.invert_yaxis()
        bins=50
        self.axHistx.hist(redshift, bins=bins, histtype='stepfilled')
        self.axHisty.hist(mag, bins=bins, histtype='stepfilled',\
                          orientation='horizontal')
        loc = plticker.MaxNLocator(nbins=4,prune='lower')
        self.axHisty.xaxis.set_major_locator(loc)
        self.axHistx.yaxis.set_major_locator(loc)
        self.axHisty.set_xlabel("N",fontsize=14)
        self.axHistx.set_ylabel("N",fontsize=14)


    def savefig(self,dir_output,name,title,tight_layout=True):

        if tight_layout:
            self.f.tight_layout(rect=[0, 0.03, 1, 0.95])
        if title != "": 
            self.f.suptitle(title,fontsize=20)
        print("Saving "+dir_output+name)
        self.f.savefig(dir_output+name,dpi=150)
        plt.close()

    ################
    # light curves #
    ################

    def plot_light_curve(self,time,signal,error,survey,band,yaxis="mag"):

        if survey in self.fmt_list.keys() : fmt = self.fmt_list[survey]
        else : fmt = {"fmt":"x","markersize":5}
        ax = self.ax_list[band]
        ax.errorbar(time,signal,yerr=error,label=band,c=self.color_list[band],\
                    **fmt)
        bottom,top = ax.get_ylim()
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
            if bottom < np.max(signal): bottom = np.max(signal)+0.1
            if top > np.min(signal): top = np.min(signal)-0.1
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

    def plot_fit_sin_curve(self,x_sin,y_sin,band):

        ax = self.ax_list[band]
        ax.plot(x_sin,y_sin,label=band,c=self.color_list[band],linestyle="--",\
                linewidth=1)


def plot_magnitude_comparison(x,y,out_dir,mjd):

    plt.scatter(x,y)
    plt.savefig(out_dir+str(mjd)+".png")
    plt.close()





