import pandas as pd
import time
import SciServer
from SciServer import Authentication, LoginPortal, Config, CasJobs, SkyQuery, SciDrive, SkyServer, Files, Jobs
import requests
from io import StringIO
from astropy import units as u
import numpy as np
import warnings
import sqlite3
from astropy.coordinates import SkyCoord
import time

class DBQuery():
    def __init__(self):
        pass
    
    def run(self):
        pass

class SDSSQuery(DBQuery):
    def __init__(self):
        from credentials.sdss import User
        self.username = User.name
        self.password = User.password
        self.__setupserver()
    
    def load_sncoor(self,sncoor,sncoor_tname='sncoor',xmatch_tname='xmatchres'):
        tables = CasJobs.getTables(context="MyDB")
        tablenames = [x['Name'] for x in tables]
        for tname in [sncoor_tname,xmatch_tname]:
            if tname in tablenames:
                self.delete_table(tname)
        self.sncoor_tname = sncoor_tname
        self.xmatch_tname = xmatch_tname
        self.sncoor = sncoor
    
    def run(self,sncoor):
        self.load_sncoor(sncoor)
        self.upload_df(self.sncoor,tablename=self.sncoor_tname)
        self.query = self._generate_querystr(self.sncoor_tname,self.xmatch_tname)
        jobid = self.submit_job(self.query)
        jobstatus = self.wait_jobdone(jobid)
        self.queryresult = self.get_query_result(self.xmatch_tname)
        
    def _generate_querystr(self,sncoortable,restable):
        query = """
                SELECT *
                INTO #snlist
                FROM MyDB.{} 
                GO

                SELECT a.oid,G.objID,G.ra,G.dec,GN.distance, PZ.z as photoz, PZ.zErr as photoz_err, S.z as specz
                into mydb.{}
                FROM #snlist a
                CROSS APPLY dbo.fGetNearbyObjEq(a.meanra, a.meandec, 0.5) GN
                JOIN Galaxy as G
                ON G.objID = GN.objID
                JOIN Photoz AS PZ
                ON PZ.objID = GN.objID
                LEFT JOIN Specobj AS S
                ON PZ.objID = S.bestObjID
                ORDER BY oid
                """.format(sncoortable,restable)
        return query
        
    def __setupserver(self):
        print("Setting up SkyServer...")
        Authentication_loginName = self.username
        Authentication_loginPassword = self.password
        token1 = Authentication.login(Authentication_loginName, Authentication_loginPassword);
        token2 = Authentication.getToken()
        token3 = Authentication.getKeystoneToken()
        token4 = Authentication.token.value
        user = Authentication.getKeystoneUserWithToken(token1)
        iden = Authentication.identArgIdentifier()
        
    def upload_df(self,df,tablename='mytable'):
        result = CasJobs.uploadPandasDataFrameToTable(dataFrame=df, tableName=tablename, context="MyDB")    
        if result:
            print("Table [{}] uploaded successfully.".format(tablename))
            
    def delete_table(self,tablename):
        sql = "DROP TABLE {}".format(tablename)
        CasJobs.executeQuery(sql=sql, context="MyDB", format="pandas")
    
    def submit_job(self,querystr,context="DR15"):
        jobId = CasJobs.submitJob(sql=querystr, context=context)
        print("Query Job is submitted. JobID={}".format(jobId))
        return jobId
        
    def wait_jobdone(self,jobId):
        jobStatus = CasJobs.getJobStatus(jobId)['Status']
        while jobStatus < 2:
            time.sleep(10)
            jobStatus = CasJobs.getJobStatus(jobId)['Status']
        if jobStatus in [2,3,4]:
            raise RuntimeError("Job is canceled or failed")
        print("Job {} is finished".format(jobId))    
            
        return jobStatus
    
    def get_query_result(self,tablename):
        query = 'select * from MyDB.{}'.format(tablename)
        df = CasJobs.getPandasDataFrameFromQuery(queryString=query, context="MyDB")
        df['distance'] = df['distance']*u.arcmin
        return df
        
        
class NEDQuery(DBQuery):
    def __init__(self,maxquery=3000):
        self.queryurl = 'https://ned.ipac.caltech.edu/cgi-bin/nnd'
        self.maxquery = int(maxquery)
    
    def load_sncoor(self,sncoor,chunk_size=200):
        self.sncoor = []
        sncoor = sncoor.reset_index(drop=True)
        if len(sncoor) > self.maxquery:
            warnings.warn("Number of coordinates is > {}. This may take too long. Only matching the first {}".format(self.maxquery,self.maxquery))
            sncoor = sncoor[0:self.maxquery]
        for i in range(0,len(sncoor),chunk_size):
            start = i
            end = np.min([i+chunk_size,len(sncoor)])
            self.sncoor.append(sncoor[start:end])
    
    def generate_url_params(self,sncoor,search_radius=30):
        uplist = ''
        for i,row in sncoor.iterrows():
            uplist += '{:.6f} {:.6f}\r\n'.format(row['meanra'],row['meandec'])

        params = {"uplist": uplist,
                  "sr_arcsec": search_radius,
                  "delimiter": "bar",
                  "NO_LINKS":1,
                  "nondb":["row_count","user_inplist","user_inp_sep"],
                  "crosid":"objname",
                  "position":["ra,dec","pretype","z","zunc"],
                 }
        self.url_params = params
        
    def run(self,sncoor,chunk_size=200):
        self.load_sncoor(sncoor,chunk_size=chunk_size)
        df = pd.DataFrame()
        for sncoor_chunk in self.sncoor:
            self.generate_url_params(sncoor_chunk)
            self.make_url_query()
            res = self.get_query_result(self.url_request_result,sncoor_chunk)
            df = pd.concat([df,res],sort=False,ignore_index=True)
#             self.queryresult = self.get_query_result(self.url_request_result)
        self.queryresult = df
    
    def make_url_query(self):
        r = requests.get(self.queryurl,params=self.url_params)
        self.url_request_result = r
    
    def get_query_result(self,r,sncoor):
        table = r.text[r.text.find('<PRE>')+len('<PRE>'):r.text.find('</PRE>')]
        data = StringIO(table)

        df = pd.read_csv(data, sep="|", header=None,skiprows=4,
                         names= ['rowno_start','location','separationArcsec','host_name','host_ra','host_dec','type','specz','specz_err','rowno_end'])
        for col in df.columns:
            if isinstance(df[col][0],str):
                df[col] = [x.strip() for x in df[col]]
        df.loc[:,'location'] = df['location'].replace('',None)
        df = df.replace({'':None})
        sncoor['location'] = ['{:.6f} {:.6f}'.format(row['meanra'],row['meandec']) for i,row in sncoor.iterrows()]
        df = df[df['type']=='G'].dropna()
        df.loc[:,'specz'] = df['specz'].astype(float)
        df.loc[:,'specz_err'] = df['specz_err'].astype(float)
        df.loc[:,'separationArcsec'] = df['separationArcsec'].astype(float)
        df = df.merge(sncoor, on='location')
        df['distance'] = df['separationArcsec']*u.arcsec.to(u.arcmin)
        df = df[['oid','host_name','distance','separationArcsec','host_ra','host_dec','specz','specz_err']]
        return df
    
    
class GLADEQuery(DBQuery):
    def __init__(self,dbfile='db/glade_v23.db'):
        self.__setup_db(dbfile=dbfile)
        self.dbfile = dbfile

    def load_sncoor(self,sncoor):
        self.sncoor = sncoor
    
    def run(self,sncoor):
        self.load_sncoor(sncoor)
        self.get_catalog_coor()
        df_coor = self.get_nearby_coor(self.catalog_coor,self.sncoor,radius=30.)
        df = self.merge_catalog_info(df_coor)
        self.queryresult = df
        
    def __setup_db(self,dbfile=None):
        print("Setting up local db [{}]".format(dbfile))
        connection = sqlite3.connect(dbfile)
        self.connection = connection
        
    def get_catalog_coor(self):
        cursor = self.connection.cursor()
        query = 'select "index",ra,dec from glade_v23'
        cursor.execute(query)
        res = cursor.fetchall()
        coor = pd.DataFrame(res,columns=['index','ra','dec'])
        self.catalog_coor = coor
    
    def get_nearby_coor(self,catalog_coor,sncoor,radius=30.):
        catalog = SkyCoord(ra=catalog_coor.ra.values*u.degree, dec=catalog_coor.dec.values*u.degree)
        c = SkyCoord(ra=sncoor.meanra*u.degree, dec=sncoor.meandec*u.degree)
        print("Cross-matching GLADE catalog using astropy module search_around_sky")
        start_time = time.time()
        idxc, idxcatalog, d2d, d3d = catalog.search_around_sky(c, radius*u.arcsec)
        df = pd.DataFrame({"idxc":idxc,"idxcatalog":idxcatalog,"separationArcsec":d2d.arcsec,"distance":d2d.arcmin})
        df = df.merge(sncoor[['oid']],left_on='idxc',right_index=True)
        print("Done. Time={} minutes".format((time.time() - start_time)/60.))
        return df
    
    def merge_catalog_info(self,df_coor,colnames=['"index"','z']):
        idxcatalog_tuple = tuple(df_coor.idxcatalog.values)
        cursor = self.connection.cursor()
        if len(idxcatalog_tuple) > 1:
            query = 'SELECT {} FROM glade_v23 where "index" in {}'.format(','.join(colnames),idxcatalog_tuple)
        else:
            query = 'SELECT {} FROM glade_v23 where "index"={}'.format(','.join(colnames),idxcatalog_tuple[0])
        cursor.execute(query)
        res = cursor.fetchall()
        cols = ['index','specz']
        if len(colnames)>2:
            cols += list(colnames[2:])
        df = pd.DataFrame(res,columns=cols)
        df = df.merge(df_coor,left_on='index',right_on='idxcatalog')
        return df
    