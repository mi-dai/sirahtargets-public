import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy import coordinates
from datetime import datetime
from astropy.time import Time
from astropy.table import Table
import os
# os.environ['SIRAHPIPE_DIR'] = '/home/mi/ztf/sirah_target_pipe'

class Broker():
    def __init__(self):
        self.queryfunc = None
    
    def make_query(self,verbose=True):
        res = self.queryfunc(self.querystr)
        self.queryresult = res
        if verbose:
            print("queryresult size:", len(self.queryresult))
    
    def load_query(self,querystr):
        self.querystr = querystr
    
    def setup(self):
        pass
    
    def gen_querystr(self,**kwargs):
        pass
    
class Alerce(Broker):
    def __init__(self):
        pass
    
    def setup(self):
        import json
        import psycopg2
        credentials_file = os.path.expandvars("$SIRAHPIPE_DIR/credentials/alercereaduser_v2.json")
        with open(credentials_file) as jsonfile:
            params = json.load(jsonfile)["params"]
        conn = psycopg2.connect(dbname=params['dbname'], user=params['user'], host=params['host'], password=params['password'])
        self.connection = conn
        self.queryfunc = lambda query: pd.read_sql_query(query, self.connection)
        
    def gen_querystr(self,mjdstart=58842,mjdend=58843,gmin=16.,gmax=20.,rmin=16.,rmax=20.,qlim=1000,realtime=True,**kwargs):
        if realtime:
            query = '''
                    select objects.oid, objects.nobs, objects.meanra, objects.meandec, objects.firstmjd, objects.lastmjd,
                    objects.first_magpsf_g, objects.first_magpsf_r, objects.last_magpsf_g, objects.last_magpsf_r, objects.classearly
                    from objects 

                    where (objects.firstmjd > {:.2f} and objects.firstmjd < {:.2f})
                    and ((objects.last_magpsf_g > {:.1f} and objects.last_magpsf_g < {:.1f})
                    or (objects.last_magpsf_r > {:.1f} and objects.last_magpsf_r < {:.1f}))
                    and (objects.classearly = 19 or objects.classearly is null)

                    limit {:d}
                    '''.format(mjdstart,mjdend,gmin,gmax,rmin,rmax,qlim)
        else:
            query = '''
                    select objects.oid, objects.nobs, objects.meanra, objects.meandec, objects.firstmjd, objects.lastmjd,
                    objects.first_magpsf_g, objects.first_magpsf_r, objects.last_magpsf_g, objects.last_magpsf_r, objects.classearly
                    from objects 

                    where (objects.firstmjd > {:.2f} and objects.firstmjd < {:.2f})
                    and (objects.classearly = 19 or objects.classearly is null)

                    limit {:d}
                    '''.format(mjdstart,mjdend,qlim)
            
        self.querystr = query
    
class Antares(Broker):
    def __init__(self):
        pass
       
    def setup(self):
        from antares_client.search import search
        self.queryfunc = search
    
class Lasair(Broker):
    """
    Credentials are stored in https://jupyter.lsst.ac.uk/
    """
    def __init__(self):
        pass
    
    def setup(self):
        import mysql.connector
        try:
            from ztf import settings
        except:
            from credentials.lasair import settings
        
        msl = mysql.connector.connect(
                user=settings.DB_USER, 
                password=settings.DB_PASS, 
                host=settings.DB_HOST, database='ztf')
        self.connection = msl
        
        self.queryfunc = lambda query: self.__queryfunc(query,self.connection)
        
    def __queryfunc(self,query,connection):
        cursor = connection.cursor()
        cursor.execute(query)
        results = cursor.fetchall()
        if len(results) == 0:
            return pd.DataFrame()
        table = Table(rows=results, 
                      names=self.qcols)
        res = table.to_pandas()
        return res
    
    def gen_querystr(self,mjdstart=58842,mjdend=58843,gmin=16.,gmax=20.,rmin=16.,rmax=20.,qlim=1000,realtime=True,use_sherlock=True,**kwargs):
        if realtime:
            self.use_sherlock = use_sherlock
            if use_sherlock:
                query = '''
                        SELECT objects.objectId, objects.ramean, objects.decmean, objects.ncand,
                        objects.jdmin - 2400000.5 AS mjdmin, objects.jdmax - 2400000.5 AS mjdmax, 
                        objects.magrmax, objects.latestrmag, objects.maggmax, objects.latestgmag,
                        sherlock_classifications.classification, sherlock_crossmatches.rank,
                        sherlock_crossmatches.catalogue_object_type, sherlock_crossmatches.catalogue_object_id,
                        sherlock_crossmatches.z, sherlock_crossmatches.photoZ, sherlock_crossmatches.photoZErr,
                        sherlock_crossmatches.catalogue_table_name, sherlock_crossmatches.separationArcsec 
                        FROM objects,sherlock_classifications,sherlock_crossmatches 
                        WHERE objects.objectId = objects.objectId 
                        AND objects.primaryId = sherlock_classifications.transient_object_id 
                        AND objects.primaryId = sherlock_crossmatches.transient_object_id 
                        AND objects.jdmin - 2400000.5 > {:.2f}
                        AND objects.jdmin - 2400000.5 < {:.2f}
                        AND ((objects.latestgmag > {:.1f} AND objects.latestgmag < {:.1f}) 
                        OR (objects.latestrmag > {:.1f} AND objects.latestrmag < {:.1f})) 
                        AND sherlock_classifications.classification NOT IN ("VS" , "AGN", "CV", "BS") 
                        AND (sherlock_crossmatches.z is not null OR sherlock_crossmatches.photoZ is not null) 
                        AND sherlock_crossmatches.catalogue_object_type = 'galaxy'
                        LIMIT {:d}
                        '''.format(mjdstart,mjdend,gmin,gmax,rmin,rmax,qlim)
            else:
                query = '''
                        SELECT objects.objectId, objects.ramean, objects.decmean, objects.ncand,
                        objects.jdmin - 2400000.5 AS mjdmin, objects.jdmax - 2400000.5 AS mjdmax, 
                        objects.magrmax, objects.latestrmag, objects.maggmax, objects.latestgmag
                        FROM objects
                        WHERE objects.jdmin - 2400000.5 > {:.2f}
                        AND objects.jdmin - 2400000.5 < {:.2f}
                        AND ((objects.latestgmag > {:.1f} AND objects.latestgmag < {:.1f}) 
                        OR (objects.latestrmag > {:.1f} AND objects.latestrmag < {:.1f}))
                        LIMIT {:d}
                        '''.format(mjdstart,mjdend,gmin,gmax,rmin,rmax,qlim)
            
        else:
            query = '''
                    SELECT objects.objectId, objects.ramean, objects.decmean, objects.ncand,
                    objects.jdmin - 2400000.5 AS mjdmin, objects.jdmax - 2400000.5 AS mjdmax, 
                    objects.magrmax, objects.latestrmag, objects.maggmax, objects.latestgmag,
                    sherlock_classifications.classification, sherlock_crossmatches.rank,
                    sherlock_crossmatches.catalogue_object_type, sherlock_crossmatches.catalogue_object_id,
                    sherlock_crossmatches.z, sherlock_crossmatches.photoZ, sherlock_crossmatches.photoZErr,
                    sherlock_crossmatches.catalogue_table_name, sherlock_crossmatches.separationArcsec 
                    FROM objects,sherlock_classifications,sherlock_crossmatches 
                    WHERE objects.objectId = objects.objectId 
                    AND objects.primaryId = sherlock_classifications.transient_object_id 
                    AND objects.primaryId = sherlock_crossmatches.transient_object_id 
                    AND objects.jdmin - 2400000.5 > {:.2f}
                    AND objects.jdmin - 2400000.5 < {:.2f}
                    AND sherlock_classifications.classification NOT IN ("VS" , "AGN", "CV", "BS") 
                    AND (sherlock_crossmatches.z is not null OR sherlock_crossmatches.photoZ is not null) 
                    AND sherlock_crossmatches.catalogue_object_type = 'galaxy'
                    LIMIT {:d}
                    '''.format(mjdstart,mjdend,qlim)
            
        self.querystr = query
        if use_sherlock:
            self.qcols = ["oid", "meanra", "meandec", "nobs" ,"firstmjd", "lastmjd", "magrmax","latestrmag",
                         "maggmax","latestgmag","classification", "xmatch_rank", "xmatch_type", "xmatch_objid",
                          "z","photoz","photoz_err","xmatch_table_name","separationArcsec"]
        else:
            self.qcols = ["oid", "meanra", "meandec", "nobs" ,"firstmjd", "lastmjd", "magrmax","latestrmag",
                         "maggmax","latestgmag"]
    
    
class TNS(Broker):
    def __init__(self):
        pass
    
    def setup(self):
        from tnsScrape import getTNS
        self.queryfunc = getTNS
    
    def gen_querystr(self,mjdstart=58842,mjdend=58843,gmin=16.,gmax=20.,rmin=16.,rmax=20.,**kwargs):
        params = {}
        params['enddate'] = Time(mjdend,format='mjd').datetime.date()
        params['discoverdays'] = mjdend - mjdstart
        params['disc_mag_min'] = np.min([gmin,rmin])
        params['disc_mag_max'] = np.max([gmax,rmax])
        params['classified'] = False
        params['skip_ztf'] = True
        params['num_page'] = 500
        params.update(kwargs)
        self.querystr = params
        
    def make_query(self,verbose=True):
        res = self.queryfunc(**self.querystr)
        self.queryresult = res
        if verbose:
            print("queryresult size:", len(self.queryresult))
        
        
        
        
        