import numpy as np
import os
import query
import plot



def show_coadd(ra,dec,band,name):

    bands = ["g","r","i"]
    image_paths = {}
    for band in bands:
        tilename,image_path = query.which_coadd_cover_it(ra,dec,band,cosmos=True)
        os.system("wget %s" % image_path)
        current_image_path = image_path.split("/")[-1]
        image_paths.update({band:current_image_path})

    plot.plot_target_color(ra,dec,image_paths,name,image_size=15)
    #plot.plot_target_coadd(ra,dec,current_image_path,band,name)
    for path in image_paths.values():
        os.system("rm %s" % path)

    return 0

if __name__ == "__main__":

    show_coadd(150.179715,2.110385,"g","J100043.1+020637")
