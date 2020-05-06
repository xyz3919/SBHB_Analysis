import math
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

        plt.rc('font', family='serif',size= 14)
        self.f,self.axes  = plt.subplots(cols, rows,**args)
        self.color_list = {"g":"g","r":"orange",\
                           "i":"brown","z":"purple"}
        if cols==4 and rows ==1:
            self.ax_list = {"g":self.axes[0],"r":self.axes[1],\
                            "i":self.axes[2],"z":self.axes[3]}
        if cols == 2 and rows == 2:
            self.ax_list = {"g":self.axes[0,0],"r":self.axes[0,1],\
                            "i":self.axes[1,0],"z":self.axes[1,1]}
            self.fmt_list = {"DES":{"fmt":"o","markersize":5,"alpha":0.3},\
                    "SDSS_corr":{"fmt":"s","markersize":5,"alpha":0.3,"markerfacecolor":"None"},\
                         "PS":{"fmt":"^","markersize":5},\
                         "ZTF":{"fmt":"x","markersize":5}
                        }

    def plot(self,x,y,band,log=False):

        self.axes.plot(x,y,label=band,c=self.color_list[band])
        if log:
            self.axes.set_xscale("log")
            self.axes.set_yscale("log")
        self.axes.legend()

    def plot_histogram(self,number_array,band,dashed=False,survey=None):

        ax = self.ax_list[band]
        ax.set_yscale("log",nonposy='clip')
        if dashed:
            kwargs = {'histtype':'step','linestyle':(0, (5,1)),'color':'grey'}
        else:
            kwargs = {'facecolor':self.color_list[band],\
                    'histtype':'stepfilled','linestyle':'solid'}

        if survey == "DES": xbound,ybound = (0,159),(0.8,500)
        elif survey == "SDSS": xbound,ybound = (0,159),(0.8,500)
        else: xbound,ybound = (0,np.percentile(number_array,99)*1.05),(0.8,500)
        ax.hist(number_array,bins=40,\
                range=xbound,\
                **kwargs)
        #if band == "i" or  band == "z":
        #    ax.set_xlabel("Number of epochs",fontsize=16)
        #if band == "g" or  band == "i":
        #    ax.set_ylabel("Number of quasars",fontsize=16)
        ax.tick_params(axis='both', which='major')
        ax.set_xlim(*xbound)
        ax.set_ylim(*ybound)
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
        self.axes.set_ylim(0,0.65)
        self.axes.tick_params(axis='both', which='major')
        self.axes.legend(prop={'size': 13},loc=2,ncol=2,columnspacing=1.2,handletextpad=0.5,labelspacing=0.5)
        self.axes.set_xlim(3300,10800)

    def plot_redshift_mag(self,redshift_raw,mag_raw):

        from mpl_toolkits.axes_grid1 import make_axes_locatable
        import matplotlib.ticker as plticker

        mag = mag_raw[((mag_raw < 50) & (mag_raw > -50))]
        redshift = redshift_raw[((mag_raw < 50) & (mag_raw > -50))]

        self.axes.scatter(redshift, mag,c="black",s=6)

        # plot periodic candidates
        """
        redshift_cand = [1.736,2.525,1.295,1.765]
        mag_cand = [22.58,20.28,20.99,20.69]
        self.axes.scatter(redshift_cand,mag_cand, marker="o", c="red",s=80,label='Candidates')
        redshift_J0252 = 1.53
        mag_J0252 = 20.6
        self.axes.scatter(redshift_J0252,mag_J0252, marker="*", c="red", s=240, label='Liao et al.' )
        """
        self.axes.scatter(1.52,20.6, marker="o", c=u'#1F77B4',edgecolors='face', s=120, label="J0252")
        self.axes.scatter(1.736,22.58, marker="s", c=u'#FF7F0E',edgecolors='face', s=100, label="J0246")
        self.axes.scatter(2.525,20.28, marker="D", c=u'#2CA02C',edgecolors='face', s=90, label="J0247")
        self.axes.scatter(1.295,20.99, marker="p", c=u'#D62728',edgecolors='face', s=150, label="J0249")
        self.axes.scatter(1.765,20.69, marker="*", c=u'#9467BD', edgecolors='face',s=240, label="J0254")



        self.axes.set_xlabel("Redshift",fontsize=16)
        self.axes.set_ylabel("r-band Magnitude",fontsize=16)
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
        self.axHisty.set_xlabel("N",fontsize=16)
        self.axHistx.set_ylabel("N",fontsize=16)

        leg = self.axes.legend(scatterpoints = 1,handletextpad=0.2, fontsize=14)

        for handle, text in zip(leg.legendHandles, leg.get_texts()):
             text.set_color(handle.get_facecolor()[0])



    def savefig(self,dir_output,name,title,xlabel=None,ylabel=None,tight_layout=True,pad=0.1,w_pad=0.0,h_pad=0.0):

        rect = [0.00,0.00,1,1]
        if title != "": 
            self.f.suptitle(title,fontsize=18)
            rect[3] = 0.93
        if xlabel is not None:
            self.f.text(0.5, 1E-3, xlabel, ha='center', va="bottom",fontsize=16)
            rect[1] = 0.05
        if ylabel is not None:
            self.f.text(1E-3, 0.5, ylabel, va='center', ha="left", rotation='vertical',fontsize=16)
            rect[0] = 0.04
        if tight_layout:
            self.f.tight_layout(rect=rect,w_pad=w_pad,h_pad=h_pad,pad=pad)

        print("Saving "+dir_output+name)
        self.f.savefig(dir_output+name,dpi=200)
        plt.close()

    ################
    # light curves #
    ################

    def plot_light_curve(self,time,signal,error,survey,band,yaxis="mag"):

        if survey in self.fmt_list.keys() : fmt = self.fmt_list[survey]
        else : fmt = {"fmt":"x","markersize":5}
        ax = self.ax_list[band]
        ax.errorbar(time,signal,yerr=error,label=band,c=self.color_list[band],\
                    mec=self.color_list[band],**fmt)
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





