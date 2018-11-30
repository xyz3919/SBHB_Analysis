import os
import numpy as np
import astropy
import query
import plot
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord


class catalog:

    def __init__(self):

        self.catalog_dir = "catalog/"
        if not os.path.isdir(self.catalog_dir): os.mkdir(self.catalog_dir)

        self.temp_catalog_dir = "catalog/temp/"
        if not os.path.isdir(self.temp_catalog_dir):
            os.mkdir(self.temp_catalog_dir)

        self.raw_catalog_dir = "raw_catalog/"
        self.milliqua = "million_quasar.txt"
        self.ozdes = "OzDES.txt"
        self.DR14Q = "DR14Q_v4_4.fits"
        self.DR7Q = "dr7qso.fit"
        self.dtype = [("ra",float),("dec", float),("flag","|S4"),\
                      ("z",float),("where","|S6")]
        self.crossmatch_radius = 2.
        self.SN_half_size = 1.1
        self.file_record = open(self.catalog_dir+"record.info","w")
        self.DES_SN_ra_dec = {"C1":(54.2743,-27.1116),
                              "C2":(54.2743,-29.0884),
                              "C3":(52.6484,-28.1000),
                              "E1":(7.8744,-43.0096),
                              "E2":(9.5000,-43.9980),
                              "S1":(42.8200,0.0000),
                              "S2":(41.1944,-0.9884),
                              "X1":(34.4757,-34.4757),
                              "X2":(35.6645,-6.4121),
                              "X3":(36.4500,-4.6000)}
        self.query = query.query()


    def load_million_catalog(self):

        print ("Loading Million quasar catalog ...")
        data = np.genfromtxt(self.raw_catalog_dir+self.milliqua,\
               delimiter=",",dtype=self.dtype)
        return data

    def load_SDSS_DR14(self):

        print ("Loading SDSS DR14 quasar catalog ...")
        hdulist = fits.open(self.raw_catalog_dir+self.DR14Q)
        data = np.array(hdulist[1].data)
        extracted_data = data[["RA","DEC","Z"]]
        extracted_data.dtype.names = "ra","dec","z"
        return extracted_data

    def load_SDSS_DR7(self):

        print ("Loading SDSS DR7 quasar catalog ...")
        hdulist = fits.open(self.raw_catalog_dir+self.DR7Q)
        data = np.array(hdulist[1].data)
        extracted_data = data[["RA","DEC","z"]]
        extracted_data.dtype.names = "ra","dec","z"
        return extracted_data

    def load_OzDES_catalog(self):

        print ("Loading OzDES quasar catalog ...")
        data = np.genfromtxt(self.raw_catalog_dir+self.ozdes,delimiter=",",\
                             dtype=[("ra",float),("dec", float),("z",float)])
        return data


    def select_quasars_in_one_region(self,data,region):

        ra,dec = self.DES_SN_ra_dec[region]
        ra_upper,ra_lower = ra+self.SN_half_size/np.cos(dec/180.*np.pi),\
                            ra-self.SN_half_size/np.cos(dec/180.*np.pi)
        dec_upper,dec_lower = dec+self.SN_half_size,dec-self.SN_half_size
        selected_data = data[(data["ra"] < ra_upper) & (data["ra"] > ra_lower) \
                           & (data["dec"] < dec_upper) & (data["dec"] > dec_lower)]
        return selected_data

    def select_quasars_in_regions(self,data,regions):

        N = 1
        for field in regions:
            selected = self.select_quasars_in_one_region(\
                       data,field)
            if N == 1 : total = selected
            else: total = np.concatenate((total,selected),axis=0)
            N = N+1
        total_unique = np.unique(total)
        total_unique.sort(order='ra')
        return total_unique

    def adding_flag(self,data,flag):

        dtype = [('ra', '<f8'), ('dec', '<f8'), ('z', '<f8')]+\
                [("flag","|S15")]
        final_list = np.array([],dtype=dtype) 
        for a in range(len(data)):
            adding_list = np.array([tuple(data[a])+(flag,)],dtype=dtype)
            final_list = np.append(final_list,adding_list)
        return final_list

    def adding_second_catalog(self,data1,data2,cat_name,spec_flag,where):

        
        self.print_and_write(self.file_record,"--- Adding "+\
                             cat_name+" catalog ---")

        c1 = SkyCoord(ra=data1["ra"]*u.degree, dec=data1["dec"]*u.degree)
        c2 = SkyCoord(ra=data2["ra"]*u.degree, dec=data2["dec"]*u.degree)
        idx, d2d, d3d = c2.match_to_catalog_sky(c1) 
        # matched
        idx_1 = idx[d2d<self.crossmatch_radius*u.arcsec]
        idx_2 = np.where(d2d<self.crossmatch_radius*u.arcsec)[0]

        self.print_and_write(self.file_record,\
                             "cross-matched:"+str(len(idx_1)))
        N = 0
        for i in range(len(idx_1)):
            if "q" in data1[idx_1[i]]["flag"]:
                if "Q" in spec_flag: 
                    data1[idx_1[i]]["flag"] = "Q"
                    N = N+1
        self.print_and_write(self.file_record,"q->Q:"+str(N))

        # not matched
        idx_not_2 = d2d>self.crossmatch_radius*u.arcsec

        adding_rows = []
        for row in data2[idx_not_2]:
            adding_rows.append((row["ra"],row["dec"],spec_flag,row["z"],where))
        adding_arrays = np.array(adding_rows,dtype=self.dtype)

        self.print_and_write(self.file_record,\
                             "adding:"+str(len(adding_arrays)))

        added_data = np.append(data1,adding_arrays)

        self.print_and_write(self.file_record,\
                             "combined:"+str(len(added_data)))

        return added_data

    def combining_two_catalogs(self,data1,data2,cat_name):

        flag_number = len(data1.dtype.names)-3+1
        dtype = [('ra', '<f8'), ('dec', '<f8'), ('z', '<f8')]+\
                [("flag_"+str(i),"|S15") for i in range(flag_number)]
        final_list = np.array([],dtype=dtype)

        self.print_and_write(self.file_record,"--- Adding "+\
                             cat_name+" catalog ---")

        c1 = SkyCoord(ra=data1["ra"]*u.degree, dec=data1["dec"]*u.degree)
        c2 = SkyCoord(ra=data2["ra"]*u.degree, dec=data2["dec"]*u.degree)
        idx, d2d, d3d = c2.match_to_catalog_sky(c1)
        # matched
        idx_1 = idx[d2d<self.crossmatch_radius*u.arcsec]
        idx_2 = np.where(d2d<self.crossmatch_radius*u.arcsec)[0]

        self.print_and_write(self.file_record,\
                             "cross-matched:"+str(len(idx_1)))
        N = 0
        for i in range(len(idx_1)):
            adding_list = np.array([tuple(data1[idx_1[i]])+\
                          (data2[idx_2[i]]["flag"],)],dtype=dtype)
            final_list = np.append(final_list,adding_list)
        # in 1 but not in 2
        mask = np.ones(len(data1), dtype=bool)
        mask[idx_1] = False
        idx_not_1 = np.where(mask)[0]

        for i in range(len(idx_not_1)):
            adding_list = np.array([tuple(data1[idx_not_1[i]])+("     ",)],\
                          dtype=dtype)
            final_list = np.append(final_list,adding_list)

        # in 2 but not in 1
        idx_not_2 = np.where(d2d>self.crossmatch_radius*u.arcsec)[0]
        self.print_and_write(self.file_record,\
                             "adding:"+str(len(idx_not_2)))

        for i in range(len(idx_not_2)):
            adding_list = np.array([tuple(data2[["ra","dec","z"]]\
                          [idx_not_2[i]])+("     ",)*(flag_number-1)+\
                          (data2["flag"][idx_not_2[i]],)],\
                          dtype=dtype)
            final_list = np.append(final_list,adding_list)

        self.print_and_write(self.file_record,\
                             "combined:"+str(len(final_list)))

        return final_list

    def adding_million_quasars(self,data1,data2):

        flag_number = len(data1.dtype.names)-3+1
        dtype = [('ra', '<f8'), ('dec', '<f8'), ('z', '<f8')]+\
                [("flag_"+str(i),"|S15") for i in range(flag_number)]
        final_list = np.array([],dtype=dtype)

        self.print_and_write(self.file_record,"--- Adding "+\
                             "Million Quasar catalog ---")

        c1 = SkyCoord(ra=data1["ra"]*u.degree, dec=data1["dec"]*u.degree)
        c2 = SkyCoord(ra=data2["ra"]*u.degree, dec=data2["dec"]*u.degree)
        idx, d2d, d3d = c2.match_to_catalog_sky(c1)
        # matched
        idx_1 = idx[d2d<self.crossmatch_radius*u.arcsec]
        idx_2 = np.where(d2d<self.crossmatch_radius*u.arcsec)[0]

        self.print_and_write(self.file_record,\
                             "cross-matched:"+str(len(idx_1)))
        N = 0
        for i in range(len(idx_1)):
            adding_list = np.array([tuple(data1[idx_1[i]])+\
                          ("MQ_"+data2[idx_2[i]]["flag"].split(" ")[0]+"_"+\
                          data2[idx_2[i]]["where"],)],dtype=dtype)
            final_list = np.append(final_list,adding_list)
        # in 1 but not in 2
        mask = np.ones(len(data1), dtype=bool)
        mask[idx_1] = False
        idx_not_1 = np.where(mask)[0]

        for i in range(len(idx_not_1)):
            adding_list = np.array([tuple(data1[idx_not_1[i]])+("     ",)],\
                          dtype=dtype)
            final_list = np.append(final_list,adding_list)

        # in 2 but not in 1
        idx_not_2 = np.where(d2d>self.crossmatch_radius*u.arcsec)[0]
        self.print_and_write(self.file_record,\
                             "adding:"+str(len(idx_not_2)))

        for i in range(len(idx_not_2)):
            adding_list = np.array([tuple(data2[["ra","dec","z"]]\
                          [idx_not_2[i]])+("     ",)*(flag_number-1)+\
                          ("MQ_"+data2["flag"][idx_not_2[i]].split(" ")[0]+"_"+\
                          data2[idx_not_2[i]]["where"],)],\
                          dtype=dtype)
            final_list = np.append(final_list,adding_list)

        self.print_and_write(self.file_record,\
                             "combined:"+str(len(final_list)))

        return final_list



    def select_spec_confirmed(self,data):

        quasar_flags = ["Q"]
        #quasar_flags = ["Q","A","B","N","K"]
        #quasar_flags = ["R","X"]
        #quasar_flags = ["q"]
        true_array = np.array([False]*len(data))
        for quasar_flag in quasar_flags:
            true_array = true_array | np.array([quasar_flag in flag \
                         for flag in data["flag"]])
        confirmed_quasars = data[true_array]
        return confirmed_quasars

    def print_and_write(self,write_file,string):

        print (string)
        write_file.write(string+"\n")


    def main(self):

        # select million quasar in SN fields
        milliqua_data = self.load_million_catalog()
        milliqua_in_SN = self.select_quasars_in_regions(milliqua_data,\
                         self.DES_SN_ra_dec.keys())

        self.print_and_write(self.file_record,"million quasars in SN:"+\
                        str(len(milliqua_in_SN)))

        np.savetxt(self.temp_catalog_dir+"milliqua_DES-SN.txt",\
                   milliqua_in_SN,fmt="%f,%f,%s,%f,%s",\
                   header="ra,dec,spec_flag,z,where")

        # adding OzDES catalog
        ozdes_data = self.load_OzDES_catalog()

        self.print_and_write(self.file_record,"OzDES quasars:"+\
                             str(len(ozdes_data)))
        
        combined = self.adding_second_catalog(milliqua_in_SN,ozdes_data,\
                                              "OzDES","Q","OzDESa")
#        np.savetxt(self.temp_catalog_dir+"milliqua+OzDES_DES-SN.txt",\
#                   combined,fmt="%f,%f,%s,%f,%s",\
#                   header="ra,dec,spec_flag,z,where")


        # cross-match with DES Y3A1 coadd catalog
        N = 1
        for row in combined:
            DES_coadd_object = self.query.get_coadd_object_spread_model(row['ra'],\
                               row['dec'])
            if DES_coadd_object is not None:
                new_row = np.array([(DES_coadd_object["ra"],\
                          DES_coadd_object["dec"],row["flag"],row["z"],\
                          row["where"],DES_coadd_object["wavg_mag_psf_r"],\
                          DES_coadd_object["spread_model_i"])],\
                          dtype=self.dtype+[("mag_psf_r",float),\
                          ("spread_model_i",float)])
                if N == 1: combined_in_DES = new_row
                else: 
                    combined_in_DES = np.concatenate((combined_in_DES,\
                                      new_row),axis=0)
                N = N+1
        self.print_and_write(self.file_record,"combined x DES Y3A1coadd:"+\
                             str(len(combined_in_DES)))

        np.savetxt(self.catalog_dir+"milliqua+OzDES_SN.txt",\
                   combined_in_DES,fmt="%f,%f,%s,%f,%s,%f,%f",\
                   header="ra,dec,spec_flag,z,where,mag_psf_r,spread_model_i")

        # extract spectroscapically confirmed quasars in S1 and S2
         
        combined_in_DES = np.genfromtxt(self.catalog_dir+\
                          "milliqua+OzDES_SN.txt",delimiter=",",\
                          dtype=self.dtype+[("mag_psf_r",float),\
                          ("spread_model_i",float)])
        quasars_S1S2 = self.select_quasars_in_regions(combined_in_DES,\
                      ["S1","S2"])
        self.print_and_write(self.file_record,"quasars in S1,S2:"+\
                             str(len(quasars_S1S2)))

        spec_quasars_S1S2 = self.select_spec_confirmed(quasars_S1S2) 
        self.print_and_write(self.file_record,"spec quasars in S1,S2:"+\
                             str(len(spec_quasars_S1S2)))
        np.savetxt(self.catalog_dir+"spec_quasars_S1S2.txt",\
                   spec_quasars_S1S2,fmt="%f,%f,%s,%f,%s,%f,%f",\
                   header="ra,dec,spec_flag,z,where,mag_psf_r,spread_model_i")

        f1 = plot.plot()
        f1.plot_ra_dec(spec_quasars_S1S2["ra"],spec_quasars_S1S2["dec"],\
                         "spec_S1S2")


if __name__ == '__main__':

    catalog = catalog()
    catalog.main()
