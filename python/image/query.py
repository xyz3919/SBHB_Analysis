import despydb.desdbi
import os,sys

def setup_dbh(section = "db-dessci"):
    # Setup desar queries here for later
    try:
        desdmfile = os.environ["des_services"]
    except KeyError:
        desdmfile = None
    dbh = despydb.desdbi.DesDbi(desdmfile,section)
    cur = dbh.cursor()
    return dbh,cur

def close(dbh,cur):

    dbh.close()

def which_coadd_cover_it(ra,dec,band,cosmos=False):

    dbh,cur = setup_dbh()
    if cosmos :
        coadd_table = "Y3A2_MISC_COADD"
        archive_table = "Y3A2_MISC_FILE_ARCHIVE_INFO"
    else : 
        coadd_table = "Y3A1_COADD"
        archive_table = "Y3A1_FILE_ARCHIVE_INFO"
    get_tile = "select tilename from %s where " % coadd_table +\
               " %s between RACMIN and RACMAX and " % ra+\
               " %s between DECCMIN and DECCMAX and band='%s'" % (dec,band)
    cur.execute(get_tile)
    tilename =  map(list,cur.fetchall())[0][0]
    get_path = "select 'https://desar2.cosmology.illinois.edu/DESFiles/desarchive/"+\
               "'||path||'/'||filename||compression"+\
               " from %s where filename in " % archive_table +\
               "(select filename from %s " % coadd_table+\
               " where %s between RACMIN and RACMAX and "% ra +\
               " %s between DECCMIN and DECCMAX and band='%s')" % (dec,band)
    print get_path
    cur.execute(get_path)
    info_path =  map(list,cur.fetchall())
    print info_path
    return tilename,info_path[0][0]


