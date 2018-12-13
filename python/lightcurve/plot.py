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
        if cols == 2 and rows == 2:
            self.ax_list = {"g":self.axes[0,0],"r":self.axes[0,1],\
                            "i":self.axes[1,0],"z":self.axes[1,1]}

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
            ax.set_xlabel("Number of epochs",size=14)
        if band == "g" or  band == "i":
            ax.set_ylabel("Number of quasars",size=14)
        ax.set_xlim(0,np.percentile(number_array,99)*1.05)
        ax.text(0.95, 0.95, band,size=14,color=self.color_list[band],\
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

    def savefig(self,dir_output,name,title):

        self.f.tight_layout(rect=[0, 0.03, 1, 0.95])
        self.f.suptitle(title,fontsize=20)
        print("Saving "+dir_output+name)
        self.f.savefig(dir_output+name,dpi=150)
        plt.close()

def plot_magnitude_comparison(x,y,out_dir,mjd):

    plt.scatter(x,y)
    plt.savefig(out_dir+str(mjd)+".png")
    plt.close()





