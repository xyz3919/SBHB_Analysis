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




