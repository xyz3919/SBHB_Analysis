import matplotlib
matplotlib.use('Agg')
import numpy as np
import os
from matplotlib import pyplot as plt


class plot:

    def __init__(self):

        self.plot_dir = "plot/"
        if not os.path.isdir(self.plot_dir): os.mkdir(self.plot_dir)

    def plot_ra_dec(self,ra,dec,name):

        plt.scatter(ra,dec,s=2)
        plt.xlabel("R.A.(deg)")
        plt.ylabel("Decl.(deg)")
        plt.savefig(self.plot_dir+"ra_dec_"+name+".png")
        plt.close()

    def plot_spread_model(self,spread_models,band,name):

        plt.hist(spread_models,bins=50,facecolor='green')
        plt.xlabel("SPREAD_MODEL_"+band,size=14)
        plt.ylabel("N",size=14)
        plt.savefig(self.plot_dir+"spread_model_i_"+name+".png")
        plt.close()

    def plot_spread_model_vs_err(self,spread_models,errs,band,name):

        plt.scatter(errs,spread_models,c="black",marker=".",s=1)
        plt.xlabel("SPREAD_MODEL_ERR_"+band,size=14)
        plt.ylabel("SPREAD_MODEL",size=14)
        plt.savefig(self.plot_dir+"spread_model_i_"+name+".png")
        plt.close()




