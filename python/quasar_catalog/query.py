import despydb.desdbi
import numpy as np
import os,sys
from astroquery.vizier import Vizier

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
                            ("spread_model_r",float),\
                            ("spread_model_err_r",float),\
                            ("mag_psf_r",float),("mag_auto_r",float)]

    def get_coadd_object_spread_model(self,ra,dec):

        dec_radian = dec*np.pi/180.
        ra_upper = (ra+2/3600./np.cos(dec_radian))
        ra_lower = (ra-2/3600./np.cos(dec_radian))
        dec_upper = dec+2/3600.
        dec_lower = dec-2/3600.
        get_object = "select RA,DEC,SPREAD_MODEL_R,SPREADERR_MODEL_R,WAVG_MAG_PSF_R,MAG_AUTO_R from Y3A1_coadd_object_summary where RA between :ra_lower and :ra_upper and DEC between :dec_lower and :dec_upper and flags_g = 0 and flags_r = 0 and flags_i = 0 and flags_z = 0 and flags_y = 0 and imaflags_iso_g = 0 and imaflags_iso_r = 0 and imaflags_iso_i = 0 and imaflags_iso_z = 0 and imaflags_iso_y = 0"
        #get_object = "select RA,DEC,SPREAD_MODEL_I,SPREADERR_MODEL_I,WAVG_MAG_PSF_I,MAG_AUTO_I from Y3A1_coadd_object_summary where RA between :ra_lower and :ra_upper and DEC between :dec_lower and :dec_upper"
        self.cur.execute(get_object,ra_lower=ra_lower,ra_upper=ra_upper,dec_lower=dec_lower,dec_upper=dec_upper)

        info_list = self.cur.fetchall()
        if len(info_list) == 1:
            return np.array(info_list,dtype=self.dtype_coadd)
        else:
            return None

    def get_OzDES_DR1_AGNs(self):

        print("... Querying the OzDES DR1 ...")
        Vizier.ROW_LIMIT = -1
        objects_OzDES = Vizier.get_catalogs(\
                        catalog="J/MNRAS/472/273/ozdesdr1")[0]
        objects_OzDES_z = objects_OzDES[((objects_OzDES["Flag"] == 3) |\
                                        (objects_OzDES["Flag"] == 4)) &\
                                        (objects_OzDES["z"] > 0.000)]
        true_array = np.array([],dtype=bool)

        print("... Getting clean AGNs from OzDES DR1  ...")
        for row in objects_OzDES_z:
            if "AGN" in row["types"]:
                if "BrightStar" not in row["types"]:
                    true_array = np.append(true_array,True)
                else: true_array = np.append(true_array,False)
            else: true_array = np.append(true_array,False)
        objects_OzDES = np.array(objects_OzDES_z[true_array][["RAJ2000",\
                        "DEJ2000","z"]])
        objects_OzDES.dtype.names = "ra","dec","z"
        return objects_OzDES 



if __name__ == "__main__":

    q = query()
    q.get_coadd_object_spread_model(40.3157,-0.749645)
