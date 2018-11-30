import os 

def degtohexname(ra,dec):

    ra = ra/15.
    rah=int(ra)
    ram=int(60*(ra-rah))
    ras=3600.*(ra-rah-(ram/60.))

    if (dec >= 0.0):
        decd=int(dec)
        decm=int(60.*(dec-decd))
        decs=3600.*(dec-decd-(decm/60.))
        return "J%02d%02d%05.2f+%02d%02d%04.1f" % (rah,ram,ras,decd,decm,decs)
    else:
        decd=-int(dec)
        dec=-dec
        decm=int(60.*(dec-decd))
        decs=3600.*(dec-decd-(decm/60.))
        dec=-dec        
        return "J%02d%02d%05.2f-%02d%02d%04.1f" % (rah,ram,ras,decd,decm,decs)


def download_file(url,path):
    
    filename = url.split("/")[-1]
    if not os.path.isfile(path):
        os.system("wget "+url)
        os.system("mv "+filename+" "+path)

def create_dir(directory):

    if not os.path.exists(directory):
        os.makedirs(directory)

def print_and_write(write_file,string):

    print (string)
    f = open(write_file,"a")
    f.write(string+"\n")
    f.close()
