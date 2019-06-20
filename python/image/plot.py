from matplotlib import pyplot as plt
import astropy
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.nddata.utils import Cutout2D
import numpy as np


def plot_target_coadd(ra,dec,filepath,band,name):


    hdu = fits.open(filepath)[1]
    wcs = WCS(hdu.header)
    rms = np.std(hdu.data)
    size = (10*u.arcsec, 10*u.arcsec)
    position = SkyCoord(ra*u.degree , dec*u.degree)
    cutout = Cutout2D(hdu.data, position, size, wcs=wcs)

    ax = plt.subplot(projection=cutout.wcs)
    plt.grid(color='white', ls='solid')
    plt.xlabel("R.A.")
    plt.ylabel("Dec.")
    plt.title(name)
    ax.imshow(cutout.data, vmin=0, vmax=4*rms )
    ax.contour(cutout.data, levels=np.logspace(np.log10(0.2*rms), np.log10(200*rms), 40), colors='white', alpha=0.5)
#    plt.savefig("figure/"name+".png")
    plt.savefig("figures/coadd_"+band+".png")
    plt.close()

def plot_target_color(ra,dec,filepaths,name,image_size=10):

    from astropy.visualization import make_lupton_rgb

    #image_size = 10 # "
    cutout_b = cutout_image(filepaths["g"],ra,dec,image_size)
    cutout_g = cutout_image(filepaths["r"],ra,dec,image_size)
    cutout_r = cutout_image(filepaths["i"],ra,dec,image_size)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rgb = make_lupton_rgb(cutout_r.data, cutout_g.data, cutout_b.data,stretch=10)
    plt.imshow(rgb, origin='lower')
    add_compass(ax,color='white')
    add_scale_bar(ax,1,image_size)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    plt.xticks([])
    plt.yticks([])
    plt.title(name+" (DES color image)")
    plt.savefig("coadd_color.png",bbox_inches='tight',dpi=100)
    plt.close()

def cutout_image(filename,ra,dec,size):

    hdu = fits.open(filename)[1]
    wcs = WCS(hdu.header)
    image_size = (size*u.arcsec, size*u.arcsec)
    position = SkyCoord(ra*u.degree , dec*u.degree)
    cutout = Cutout2D(hdu.data, position, image_size, wcs=wcs)

    return cutout


def add_compass(ax,color='black'):

    compass_orig = ax.get_xlim()[1]-(ax.get_xlim()[1]-ax.get_xlim()[0])*0.05,ax.get_ylim()[0]+(ax.get_ylim()[1]-ax.get_ylim()[0])*0.05
    compass_length = (ax.get_ylim()[1]-ax.get_ylim()[0])*0.1
    ax.arrow(compass_orig[0],compass_orig[1],-compass_length,0, color=color,head_width=0.1*compass_length,head_length=0.15*compass_length)
    ax.arrow(compass_orig[0],compass_orig[1],0,compass_length, color=color,head_width=0.1*compass_length,head_length=0.15*compass_length)
    ax.text(compass_orig[0]-compass_length*1.35,compass_orig[1],"E",verticalalignment='center', horizontalalignment='right',color=color)
    ax.text(compass_orig[0],compass_orig[1]+compass_length*1.35,"N",verticalalignment='bottom', horizontalalignment='center',color=color)

def add_scale_bar(ax,scale_bar_size,image_size):

    scale_bar_orig = ax.get_xlim()[0]+(ax.get_xlim()[1]-ax.get_xlim()[0])*0.05,ax.get_ylim()[0]+(ax.get_ylim()[1]-ax.get_ylim()[0])*0.05
    scale_bar_length = float(scale_bar_size)/image_size*(ax.get_ylim()[1]-ax.get_ylim()[0])
    ax.arrow(scale_bar_orig[0],scale_bar_orig[1],scale_bar_length,0, head_length=scale_bar_length*0.01, color='white')
    ax.text(scale_bar_orig[0]+scale_bar_length*0.5,scale_bar_orig[1],str(scale_bar_size)+'"',verticalalignment='bottom', horizontalalignment='center',color='white')

