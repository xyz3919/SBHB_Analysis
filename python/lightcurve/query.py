import despydb.desdbi
import numpy as np
import os,sys

import string,urllib,copy
from math import log10, radians, pi,cos,sin
import time, traceback, datetime
import threading
import StringIO

#############
# query DES #
#############

class query_DES:

    def __init__(self):

        self.dtype_single = [("ra",float),("dec",float),("mjd_obs",float),\
                            ("flux_psf",float),("flux_err_psf",float),\
                            ("flux_auto",float),("flux_err_auto",float),\
                            ("band","|S1")]
        self.dtype_coadd = [("ra",float),("dec",float),("mag_psf_g",float),\
                            ("mag_psf_r",float),("mag_psf_i",float),\
                            ("mag_psf_z",float),("mag_psf_Y",float)]
        try:
            desdmfile = os.environ["des_services"]
        except KeyError:
            desdmfile = None
        self.dbh = despydb.desdbi.DesDbi(desdmfile,"db-dessci")
        self.cur = self.dbh.cursor()

        self.dbh_oper = despydb.desdbi.DesDbi(desdmfile,"db-desoper")
        self.cur_oper = self.dbh_oper.cursor()

    def get_nearby_single_epoch_objects(self,ra,dec,radius):

        dec_radian = dec*np.pi/180.
        ra_upper = (ra+radius/3600./np.cos(dec_radian))
        ra_lower = (ra-radius/3600./np.cos(dec_radian))
        dec_upper = dec+radius/3600.
        dec_lower = dec-radius/3600.
        # get Y1-Y5 objects (with y4a1_v1.5 zeropoint calibration)
        get_list = "with SNVAR_TEMP_1 as (select CATALOGNAME, MAG_ZERO, SIGMA_MAG_ZERO, MAG_ONE, SIGMA_MAG_ONE, MJD_OBS from  Y6A1_EXPOSURE join Y6A1_ZEROPOINT on Y6A1_EXPOSURE.expnum = Y6A1_ZEROPOINT.expnum where  Y6A1_ZEROPOINT.flag < 16 and Y6A1_ZEROPOINT.source = 'FGCM' and Y6A1_ZEROPOINT.version =:version and MJD_OBS>56400 and EXPTIME > 30 order by CATALOGNAME ) select RA, DEC, MJD_OBS, Y6A1_FINALCUT_OBJECT.FLUX_PSF*POWER(10,-0.4*MAG_ZERO+9) as FLUX_PSF, SQRT(POWER(1.09*SIGMA_MAG_ZERO*Y6A1_FINALCUT_OBJECT.FLUX_PSF, 2) + POWER(Y6A1_FINALCUT_OBJECT.FLUXERR_PSF, 2))*POWER(10,-0.4*MAG_ZERO+9) as FLUXERR_PSF, Y6A1_FINALCUT_OBJECT.FLUX_AUTO*POWER(10,-0.4*MAG_ZERO+9) as FLUX_AUTO, SQRT(POWER(1.09*SIGMA_MAG_ZERO*Y6A1_FINALCUT_OBJECT.FLUX_AUTO, 2) + POWER(Y6A1_FINALCUT_OBJECT.FLUXERR_AUTO, 2))*POWER(10,-0.4*MAG_ZERO+9) as FLUXERR_AUTO, BAND from  Y6A1_FINALCUT_OBJECT, SNVAR_TEMP_1 where CATALOGNAME = FILENAME and Flags = 0 and imaflags_iso = 0 and RA between :ra_lower and :ra_upper and DEC between :dec_lower and :dec_upper  order by RA, DEC"
        self.cur.execute(get_list,version='y5a1_v1.1',ra_lower=ra_lower,ra_upper=ra_upper,dec_lower=dec_lower,dec_upper=dec_upper)
        info_list1 = self.cur.fetchall()
        # get SV objects (with v2.0 zeropoint calibration)
        get_list = "with SNVAR_TEMP_1 as (select CATALOGNAME, MAG_ZERO, SIGMA_MAG_ZERO, MAG_ONE, SIGMA_MAG_ONE, MJD_OBS from  Y4A1_EXPOSURE join Y4A1_ZEROPOINT on Y4A1_EXPOSURE.expnum = Y4A1_ZEROPOINT.expnum where  Y4A1_ZEROPOINT.flag < 16 and Y4A1_ZEROPOINT.source = 'FGCM' and Y4A1_ZEROPOINT.version = :version and MJD_OBS<56400 and EXPTIME > 30 order by CATALOGNAME ) select RA, DEC, MJD_OBS, Y4A1_FINALCUT_OBJECT.FLUX_PSF*POWER(10,-0.4*MAG_ZERO+9) as FLUX_PSF, SQRT(POWER(1.09*SIGMA_MAG_ZERO*Y4A1_FINALCUT_OBJECT.FLUX_PSF, 2) + POWER(Y4A1_FINALCUT_OBJECT.FLUXERR_PSF, 2))*POWER(10,-0.4*MAG_ZERO+9) as FLUXERR_PSF, Y4A1_FINALCUT_OBJECT.FLUX_AUTO*POWER(10,-0.4*MAG_ZERO+9) as FLUX_AUTO, SQRT(POWER(1.09*SIGMA_MAG_ZERO*Y4A1_FINALCUT_OBJECT.FLUX_AUTO, 2) + POWER(Y4A1_FINALCUT_OBJECT.FLUXERR_AUTO, 2))*POWER(10,-0.4*MAG_ZERO+9) as FLUXERR_AUTO, BAND from  Y4A1_FINALCUT_OBJECT, SNVAR_TEMP_1 where CATALOGNAME = FILENAME and Flags = 0 and imaflags_iso = 0 and RA between :ra_lower and :ra_upper and DEC between :dec_lower and :dec_upper  order by RA, DEC"
        self.cur.execute(get_list,version="v2.0",ra_lower=ra_lower,ra_upper=ra_upper,dec_lower=dec_lower,dec_upper=dec_upper)
        info_list2 = self.cur.fetchall()
        info_list = info_list1+info_list2
        if len(info_list) > 0:
            objects = np.array(info_list,dtype=self.dtype_single)
            return objects
        else:
            return None

    def get_nearby_coadd_objects(self,ra,dec,radius,addition=""):

        dec_radian = dec*np.pi/180.
        ra_upper = (ra+radius/3600./np.cos(dec_radian))
        ra_lower = (ra-radius/3600./np.cos(dec_radian))
        dec_upper = dec+radius/3600.
        dec_lower = dec-radius/3600.
        get_list = "select RA,DEC,WAVG_MAG_PSF_G,WAVG_MAG_PSF_R,WAVG_MAG_PSF_I,WAVG_MAG_PSF_Z,WAVG_MAG_PSF_Y from Y3A2_COADD_OBJECT_SUMMARY where WAVG_MAG_PSF_I between 15 and 28  and RA between :ra_lower and :ra_upper and DEC between :dec_lower and :dec_upper "+addition
        self.cur.execute(get_list,ra_lower=ra_lower,ra_upper=ra_upper,dec_lower=dec_lower,dec_upper=dec_upper)
        info_list = self.cur.fetchall()
        if len(info_list) > 0:
            objects = np.array(info_list,dtype=self.dtype_coadd)
            return objects
        else:
            return None

    def get_filename_list_from_tag(self,ra,dec,tag):

        dec_radian = dec*np.pi/180.
        ra_upper = (ra+10.0/3600./np.cos(dec_radian))
        ra_lower = (ra-10.0/3600./np.cos(dec_radian))
        dec_upper = dec+10.0/3600.
        dec_lower = dec-10.0/3600.
        get_list = "select distinct o.filename,o.expnum,f.id from SE_object o ,proctag p, pfw_attempt f where RA between :ra_lower and :ra_upper and DEC between :dec_lower and :dec_upper and 'D00'||to_char(o.expnum)=f.unitname and f.id=p.pfw_attempt_id and p.tag=:tag and o.band not in ('u','Y') order by o.filename"
        self.cur_oper.execute(get_list,ra_lower=ra_lower,ra_upper=ra_upper,dec_lower=dec_lower,dec_upper=dec_upper,tag=tag)
        info_list =  self.cur_oper.fetchall()
        self.dtype_list_filename = [("filename","|S41"),("expnum",int),\
                                    ("id",int)]
        if len(info_list) > 0:
            filename_list = np.array(info_list,dtype=self.dtype_list_filename)
            return filename_list
        else:
            return None

    def get_firstcut_objects_from_filename(self,filename,tag,star=False):

        if star:
            get_list = "select distinct ra,dec,MJD_OBS,FLUX_PSF,FLUXERR_PSF,FLUX_AUTO,FLUXERR_AUTO,o.BAND from SE_object o ,proctag p, pfw_attempt f, exposure e where e.expnum=o.expnum and 'D00'||to_char(o.expnum)=f.unitname and f.id=p.pfw_attempt_id and p.tag='"+str(tag)+"' and o.filename=:filename and Flags = 0 and imaflags_iso = 0 and spread_model < 0.005"
        else:
            get_list = "select distinct ra,dec,MJD_OBS,FLUX_PSF,FLUXERR_PSF,FLUX_AUTO,FLUXERR_AUTO,o.BAND from SE_object o ,proctag p, pfw_attempt f, exposure e where e.expnum=o.expnum and 'D00'||to_char(o.expnum)=f.unitname and f.id=p.pfw_attempt_id and p.tag='"+str(tag)+"' and o.filename=:filename and Flags = 0 and imaflags_iso = 0"
        self.cur_oper.execute(get_list,filename=filename)
        info_list = self.cur_oper.fetchall()
        objects = np.array(info_list,dtype=self.dtype_single)
        if len(objects) == 0 : print("Find 0 objects")
        return objects

#######################
# query SDSS Strip 82 #
#######################

global limit_per_min, recent_query_list, do_limit
limit_per_min = 19.99
recent_query_list = []
do_limit = True

class Stripe82:
        """
        Written by J. Bloom (jsbloom@astro.berkeley.edu), Aug. 2010.
        Version 1.0

        Edited by Yu-Ching(Tony) Chen (ycchen@illinois.ede), Oct 2018.
        """

        dr_url="http://cas.sdss.org/stripe82/en/tools/search/x_sql.asp"
        formats = ['csv','xml','html']
        def_fmt = "csv"
        verbose = True

        def __init__(self):
                self.dtype = [("ra",float),("dec",float),\
                              ("mag_psf_g",float),("mag_err_psf_g",float),\
                              ("mag_psf_r",float),("mag_err_psf_r",float),\
                              ("mag_psf_i",float),("mag_err_psf_i",float),\
                              ("mag_psf_z",float),("mag_err_psf_z",float),\
                              ("mjd_g",float),\
                              ("mjd_r",float),("mjd_i",float),("mjd_z",float)]

        def _filtercomment(self,sql):
                "Get rid of comments starting with --"
                fsql = ''
                for line in sql.split('\n'):
                        fsql += line.split('--')[0] + ' ' + os.linesep;
                return fsql

        def q(self,ra,dec,dist=4.0,clobber=False,idn=0):

                sss = """SELECT p.ra,p.dec,
                                                p.psfmag_g,p.psfMagErr_g,
                                                p.psfmag_r,p.psfMagErr_r,
                                                p.psfmag_i,p.psfMagErr_i,
                                                p.psfmag_z,p.psfMagErr_z,
                                f.mjd_r,f.mjd_g,f.mjd_i,f.mjd_z
                             FROM PhotoObjAll p
                                 JOIN field as f on p.fieldid = f.fieldid
                         WHERE ra between %f and %f
                           and dec between %f and %f
                           and p.parentID = 0
                           order by f.mjd_r
                      """ % ( ra - cos(radians(dec))*0.5*dist/3600.0, ra + cos(radians(dec))*0.5*dist/3600.0,\
               dec - 0.5*dist/3600.0, dec + 0.5*dist/3600.0)
                fff = self._query(sss)
                line = fff.readline()
                if line.startswith("ERROR") or line.startswith("No objects"):
                        if self.verbose:
                                print line
                                print fff.readlines()
                        return
                ## get the keys
                kkk = line.split(",")
                print kkk
                """
                f = open(outname,"w")
                f.write("# stripe82 phot for (ra,dec) = (%f,%f)\n" % (ra,dec))
                f.write("# time run %s (utc)\n" % str(datetime.datetime.utcnow()))
                f.write(",".join(kkk))
                """
                ooo = fff.readlines()
                list_objects = [ map(float,o[:-1].split(",")) for o in ooo]
                tuple_objects = map(tuple, list_objects)
                objects =  np.array(tuple_objects,dtype=self.dtype)

                print " .... %i data points" % len(objects)
                #f.writelines(ooo)
                return objects

        def q_orig(self,ra,dec,dist=5.0,clobber=False,idn=0):
                outname = self.outdir + "sdss_%i_%f%s%f.phot" % (idn,ra,"-" if dec < 0 else "+",abs(dec))
                if not clobber and os.path.exists(outname):
                        print "file %s already exists" % outname
                        return

                sss = """SELECT p.objid,p.ra,p.dec,dbo.fdistancearcmineq(p.ra,p.dec,%f,%f)*60 as darcsec,
                                p.dered_r,
                                                p.psfmag_u,p.psfMagErr_u,p.dered_u,p.u,p.err_u,
                                                p.psfmag_g,p.psfMagErr_g,p.dered_g,p.g,p.err_g,
                                                p.psfmag_r,p.psfMagErr_r,p.dered_r,p.r,p.err_r,
                                                p.psfmag_i,p.psfMagErr_i,p.dered_i,p.i,p.err_i,
                                p.psfmag_z,p.psfMagErr_z,p.dered_z,p.z,p.err_z,
                                                p.petroR90_r,p.petroR90Err_r,
                                case
                                                  when p.parentID != 0 then "True"
                                                  else "False"
                                                end as "deblended",
                                f.mjd_r,f.mjd_g,f.mjd_i,f.nobjects,f.nstars
                             FROM PhotoObjAll p
                                 JOIN field as f on p.fieldid = f.fieldid
                         WHERE ra between %f and %f
                           and dec between %f and %f
                           order by f.mjd_r
                      """ % (ra, dec, ra - cos(radians(dec))*0.5*dist/3600.0, ra + cos(radians(dec))*0.5*dist/3600.0,\
               dec - 0.5*dist/3600.0, dec + 0.5*dist/3600.0)
                fff = self._query(sss)
                line = fff.readline()
                if line.startswith("ERROR") or line.startswith("No objects"):
                        if self.verbose:
                                print line
                                print fff.readlines()
                        return
                ## get the keys
                kkk = line.split(",")
                print kkk
                f = open(outname,"w")
                f.write("# stripe82 phot for (ra,dec) = (%f,%f)\n" % (ra,dec))
                f.write("# time run %s (utc)\n" % str(datetime.datetime.utcnow()))
                f.write("#" + ",".join(kkk))
                ooo = fff.readlines()
                print " .... %i data points" % len(ooo)
                f.writelines(ooo)
                return

        def runnat(self,f="qso_not_in_sesar.txt"):
                a= numpy.loadtxt(f)
                for i in range(len(a[:,0])):
                        ra,dec = a[i,0],a[i,1]
                        print i, ra, dec
                        self.q(ra,dec)

        def _query(self,sql,url=dr_url,fmt=def_fmt,wait_period=62):
                "Run query and return file object"
                global limit_per_min, recent_query_list, do_limit
                if do_limit:
                        recent_query_list.append((time.time(),sql))
                        ## cull out all the old calls
                        recent_query_list = [x for x in recent_query_list if time.time() - x[0] < wait_period]
                        if len(recent_query_list) > limit_per_min:
                                ## ug, we've got to wait
                                tmp = [time.time() - x[0] for x in recent_query_list if time.time() - x[0] < wait_period]
                                wait = wait_period - max(tmp) + 1
                                if self.verbose:
                                        print "Date: %s Query length is %i in the last %f sec" % (str(datetime.datetime.now()),len(recent_query_list) - 1 , wait_period)
                                        print "waiting %f sec, so as not to block the SDSS query %s" % (wait,sql)
                                time.sleep(wait)

                fsql = self._filtercomment(sql)
                params = urllib.urlencode({'cmd': fsql, 'format': fmt})
                try:
                        return urllib.urlopen(url+'?%s' % params)
                except:
                        print "TRIED: " + url+'?%s' % params
                        print "EXCEPT: sdss.py._query()"
                        return StringIO.StringIO() # This is an empty filehandler

###############
# query SDSS  #
###############

class query_SDSS:

    """>> sqlcl << command line query tool by Tamas Budavari <budavari@jhu.edu>
    Usage: sqlcl [options] sqlfile(s)

    Options:
        -s url     : URL with the ASP interface (default: pha)
        -f fmt     : set output format (html,xml,csv - default: csv)
        -q query   : specify query on the command line
        -l         : skip first line of output with column names
        -v         : verbose mode dumps settings in header
        -h         : show this message"""


    def __init__(self):

        self.formats = ['csv','xml','html']
        public_url='http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'
        self.default_url=public_url
        self.default_fmt='csv'


    def usage(self,status, msg=''):
        "Error message and usage"
        print __doc__
        if msg:
            print '-- ERROR: %s' % msg
        sys.exit(status)

    def filtercomment(self,sql):
        "Get rid of comments starting with --"
        import os
        fsql = ''
        for line in sql.split('\n'):
            fsql += line.split('--')[0] + ' ' + os.linesep;
        return fsql

    def query(self,sql,url,fmt):
        "Run query and return file object"
        import urllib
        fsql = self.filtercomment(sql)
        params = urllib.urlencode({'cmd': fsql, 'format': fmt})
        return urllib.urlopen(url+'?%s' % params)

    def write_header(self,ofp,pre,url,qry):
        import  time
        ofp.write('%s SOURCE: %s\n' % (pre,url))
        ofp.write('%s TIME: %s\n' % (pre,time.asctime()))
        ofp.write('%s QUERY:\n' % pre)
        for l in qry.split('\n'):
            ofp.write('%s   %s\n' % (pre,l))

    def main(self,argv):
        "Parse command line and do it..."
        import os, getopt, string

        queries = []
        url = os.getenv("SQLCLURL",self.default_url)
        fmt = self.default_fmt
        writefirst = 1
        verbose = 0

        # Parse command line
        try:
            optlist, args = getopt.getopt(argv[0:],'s:f:q:vlh?')
        except getopt.error, e:
            self.usage(1,e)

        for o,a in optlist:
            if   o=='-s': url = a
            elif o=='-f': fmt = a
            elif o=='-q': queries.append(a)
            elif o=='-l': writefirst = 0
            elif o=='-v': verbose += 1
            else: self.usage(0)

        if fmt not in self.formats:
            self.usage(1,'Wrong format!')

        # Enqueue queries in files
        for fname in args:
            try:
                queries.append(open(fname).read())
            except IOError, e:
                self.usage(1,e)

        # Run all queries sequentially
        for qry in queries:
            ofp = sys.stdout
            if verbose:
                self.write_header(ofp,'#',url,qry)
            file = self.query(qry,url,fmt)
            # Output line by line (in case it's big)
            line = file.readline()
            if line.startswith("ERROR"): # SQL Statement Error -> stderr
                ofp = sys.stderr
                return None,None
            '''
            if writefirst:
                ofp.write(string.rstrip(line)+os.linesep)
            line = file.readline()
            while line:
                ofp.write(string.rstrip(line)+os.linesep)
                line = file.readline()
            '''
            column = file.readline()
            data = file.read()
            print column,data
            if not data: return column,None
            return column,data

####################
# query Pan-STARRS #
####################

# NOT YET FINISHED !!!

class query_panstarrs():

    """ 
    from the website 'Query Pan-STARRS DR2 catalog using MAST API'
    url:http://ps1images.stsci.edu/ps1_dr2_api.html
    """

    def ps1cone(ra,dec,radius,table="mean",release="dr1",\
                format="csv",columns=None,\
                baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs",\
                verbose=False,**kw):
        """Do a cone search of the PS1 catalog
    
        Parameters
        ----------
        ra (float): (degrees) J2000 Right Ascension
        dec (float): (degrees) J2000 Declination
        radius (float): (degrees) Search radius (<= 0.5 degrees)
        table (string): mean, stack, or detection
        release (string): dr1 or dr2
        format: csv, votable, json
        columns: list of column names to include (None means use defaults)
        baseurl: base URL for the request
        verbose: print info about request
        **kw: other parameters (e.g., 'nDetections.min':2)
        """
    
        data = kw.copy()
        data['ra'] = ra
        data['dec'] = dec
        data['radius'] = radius
        return ps1search(table=table,release=release,format=format,\
                         columns=columns,baseurl=baseurl, verbose=verbose,\
                         **data)

    def ps1search(table="mean",release="dr1",format="csv",columns=None,
                  baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs",\
                  verbose=False, **kw):
        """Do a general search of the PS1 catalog (possibly without ra/dec/radius)
        Parameters
        ----------
        table (string): mean, stack, or detection
        release (string): dr1 or dr2
        format: csv, votable, json
        columns: list of column names to include (None means use defaults)
        baseurl: base URL for the request
        verbose: print info about request
        **kw: other parameters (e.g., 'nDetections.min':2).  Note this is required!
        """

        data = kw.copy()
        if not data:
            raise ValueError("You must specify some parameters for search")
        self.checklegal(table,release)
        if format not in ("csv","votable","json"):
            raise ValueError("Bad value for format")
        url = "{baseurl}/{release}/{table}.{format}".format(**locals())
        if columns:
            # check that column values are legal
            # create a dictionary to speed this up
            dcols = {}
            for col in ps1metadata(table,release)['name']:
                dcols[col.lower()] = 1
            badcols = []
            for col in columns:
                if col.lower().strip() not in dcols:
                    badcols.append(col)
            if badcols:
                raise ValueError('Some columns not found in table: {}'.format(', '.join(badcols)))
            # two different ways to specify a list of column values in the API
            # data['columns'] = columns
            data['columns'] = '[{}]'.format(','.join(columns))

    # either get or post works
    #    r = requests.post(url, data=data)
        r = requests.get(url, params=data)
    
        if verbose:
            print(r.url)
        r.raise_for_status()
        if format == "json":
            return r.json()
        else:
            return r.text

    def checklegal(table,release):
        """Checks if this combination of table and release is acceptable
    
        Raises a VelueError exception if there is problem
        """
    
        releaselist = ("dr1", "dr2")
        if release not in ("dr1","dr2"):
            raise ValueError("Bad value for release (must be one of {})".format(', '.join(releaselist)))
        if release=="dr1":
            tablelist = ("mean", "stack")
        else:
            tablelist = ("mean", "stack", "detection")
        if table not in tablelist:
            raise ValueError("Bad value for table (for {} must be one of {})".format(release, ", ".join(tablelist)))

    def ps1metadata(table="mean",release="dr1",\
                  baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs"):
        """Return metadata for the specified catalog and table
    
        Parameters
        ----------
        table (string): mean, stack, or detection
        release (string): dr1 or dr2
        baseurl: base URL for the request
    
        Returns an astropy table with columns name, type, description
        """
    
        checklegal(table,release)
        url = "{baseurl}/{release}/{table}/metadata".format(**locals())
        r = requests.get(url)
        r.raise_for_status()
        v = r.json()
        # convert to astropy table
        tab = Table(rows=[(x['name'],x['type'],x['description']) for x in v],
                    names=('name','type','description'))
        return tab

    def mastQuery(request):
        """Perform a MAST query.

        Parameters
        ----------
        request (dictionary): The MAST request json object
    
        Returns head,content where head is the response HTTP headers, and content is the returned data
        """
        
        server='mast.stsci.edu'
    
        # Grab Python Version 
        version = ".".join(map(str, sys.version_info[:3]))
    
        # Create Http Header Variables
        headers = {"Content-type": "application/x-www-form-urlencoded",
                   "Accept": "text/plain",
                   "User-agent":"python-requests/"+version}

        # Encoding the request as a json string
        requestString = json.dumps(request)
        requestString = urlencode(requestString)
    
        # opening the https connection
        conn = httplib.HTTPSConnection(server)

        # Making the query
        conn.request("POST", "/api/v0/invoke", "request="+requestString, headers)

        # Getting the response
        resp = conn.getresponse()
        head = resp.getheaders()
        content = resp.read().decode('utf-8')
    
        # Close the https connection
        conn.close()

        return head,content

    def resolve(name):
        """Get the RA and Dec for an object using the MAST name resolver
    
        Parameters
        ----------
        name (str): Name of object

        Returns RA, Dec tuple with position"""

        resolverRequest = {'service':'Mast.Name.Lookup',
                           'params':{'input':name,
                                     'format':'json'
                                    },
                          }
        headers,resolvedObjectString = mastQuery(resolverRequest)
        resolvedObject = json.loads(resolvedObjectString)
        # The resolver returns a variety of information about the resolved object, 
        # however for our purposes all we need are the RA and Dec
        try:
            objRa = resolvedObject['resolvedCoordinate'][0]['ra']
            objDec = resolvedObject['resolvedCoordinate'][0]['decl']
        except IndexError as e:
            raise ValueError("Unknown object '{}'".format(name))
        return (objRa, objDec)

