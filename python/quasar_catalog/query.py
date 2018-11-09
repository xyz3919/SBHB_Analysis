import despydb.desdbi
import numpy as np
import os,sys

class query:

    def __init__(self):

        try:
            desdmfile = os.environ["des_services"]
        except KeyError:
            desdmfile = None
        self.dbh = despydb.desdbi.DesDbi(desdmfile,"db-dessci")
        self.cur = self.dbh.cursor()

        self.dbh_oper = despydb.desdbi.DesDbi(desdmfile,"db-desoper")
        self.cur_oper = self.dbh_oper.cursor()

        self.dtype_coadd = [("ra",float),("dec",float),\
                            ("spread_model_r",float),("flux_psf_r",float),\
                            ("flux_auto_r",float)]

    def get_coadd_object_spread_model(self,ra,dec):

        dec_radian = dec*np.pi/180.
        ra_upper = (ra+2/3600./np.cos(dec_radian))
        ra_lower = (ra-2/3600./np.cos(dec_radian))
        dec_upper = dec+2/3600.
        dec_lower = dec-2/3600.
        get_object = "select RA,DEC,SPREAD_MODEL_R,WAVG_FLUX_PSF_R,FLUX_AUTO_R from Y3A1_coadd_object_summary where RA between :ra_lower and :ra_upper and DEC between :dec_lower and :dec_upper and flags_g = 0 and flags_r = 0 and flags_i = 0 and flags_z = 0 and flags_y = 0 and imaflags_iso_g = 0 and imaflags_iso_r = 0 and imaflags_iso_i = 0 and imaflags_iso_z = 0 and imaflags_iso_y = 0"
        self.cur.execute(get_object,ra_lower=ra_lower,ra_upper=ra_upper,dec_lower=dec_lower,dec_upper=dec_upper)

        info_list = self.cur.fetchall()
        if len(info_list) == 1:
            return np.array(info_list,dtype=self.dtype_coadd)
        else:
            return None

if __name__ == "__main__":

    q = query()
    q.get_coadd_object_spread_model(40.3157,-0.749645)
