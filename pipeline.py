import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import inspect

class SIRAHPipe():
    def __init__(self,brokers=['alerce','lasair','tns'],xmatch_catalogues=['sdss','ned','glade'],
                 selection_cuts=['zcut','ztflimcut','hostdistcut','olddetectioncut','magzcut','rbcut']):
        """
        Parameters
        ----------
        brokers
        xmatch_catalogues
        selection_cuts
        
        """
        self.Brokers = Brokers(brokers=brokers)
#         if np.any([x in brokers for x in ['alerce','tns']]):
        self.CrossMatch = CrossMatch(catalogues=xmatch_catalogues)
        self.MakeCuts = MakeCuts()
        self.selection_cuts = selection_cuts
        print("Brokers to query:",brokers)
        print("Crossmatch catalogues:",xmatch_catalogues)
        print("Cuts to apply:",selection_cuts)
            
    def run(self,realtime=True,mjdstart=58842,mjdend=58843,**kwargs):
        
        self.Brokers.run(realtime=realtime,mjdstart=mjdstart,mjdend=mjdend,**kwargs)
        self.results = pd.DataFrame()
        if 'alerce' in self.Brokers.names:
            idx = self.Brokers.names.index('alerce')
            sncoor = self.Brokers.queryresults[idx][['oid','meanra','meandec']]
            if len(sncoor) == 0:
                print("query result is empty")
            else:                              
                #crossmatch with sdss/glade first
                dblist = [x for x in ['sdss','glade'] if x in self.CrossMatch.db_dict.keys()]
                self.CrossMatch.run(sncoor=sncoor,select_db=dblist)        
                xmatch_res = self.CrossMatch.result.sort_values('distance')

                if 'ned' in self.CrossMatch.db_dict.keys():
                    #crossmatch ned for photoz in range
                    search_idx = (xmatch_res['photoz'] > 0.01) & (xmatch_res['photoz'] < 0.08) & (xmatch_res['specz'].isna())
                    if np.sum(search_idx) > 0:
                        print("Cross match NED for 0.01 < sdss_photoz < 0.08, num = ",np.sum(search_idx))
                        search_oid = xmatch_res.loc[search_idx,'oid'].drop_duplicates()
                        coor = sncoor.set_index('oid').loc[search_oid].reset_index()
                        self.CrossMatch.run(sncoor=coor,select_db=['ned'])
                        xmatch_res = pd.concat([xmatch_res,self.CrossMatch.result.sort_values('distance')],sort=False,ignore_index=True)

                    #crossmatch ned for dec<0
                    south_idx = sncoor['meandec'] < 0
                    if np.sum(south_idx) > 0:
                        print("Cross match NED for dec < 0, num = ",np.sum(south_idx))
                        coor = sncoor.loc[south_idx].drop_duplicates('oid')
                        self.CrossMatch.run(sncoor=coor,select_db=['ned'])
                        xmatch_res = pd.concat([xmatch_res,self.CrossMatch.result.sort_values('distance')],sort=False,ignore_index=True)

                xmatch_res = xmatch_res.sort_values(['oid','distance'])   
                res = self.Brokers.queryresults[idx].merge(xmatch_res,how='inner',on='oid')
                res = self._assign_redshift_from_sdss(res)
                res = self._rename_columns(res,orignames=['last_magpsf_g','last_magpsf_r'],
                                           newnames=['gmaglatest','rmaglatest'])
                res['Broker'] = 'alerce'
                
                self.results = pd.concat([self.results,res],sort=False,join='outer',ignore_index=True)
                

        if 'tns' in self.Brokers.names:
            idx = self.Brokers.names.index('tns')
            sncoor = self.Brokers.queryresults[idx][['oid','meanra','meandec']]
            if len(sncoor) == 0:
                print("query result is empty")
            else:
                dblist = [x for x in ['sdss','ned','glade'] if x in self.CrossMatch.db_dict.keys()]
                self.CrossMatch.run(sncoor=sncoor,select_db=dblist)   
                res = self.Brokers.queryresults[idx].merge(self.CrossMatch.result.sort_values('distance'),how='inner',on='oid')
                res = self._assign_redshift_from_sdss(res)
                res = self._rename_columns(res,orignames=['discoverymag'],
                                           newnames=['gmaglatest'])
                res['Broker'] = 'tns'
                res['nobs'] = -99
                res['gmaglatest'] = res['gmaglatest'].astype(float)
                res['rmaglatest'] = np.nan
                res['lastmjd'] = Time(list(res['discoverydate'].values),format='iso').mjd
                self.results = pd.concat([self.results,res],sort=False,join='outer',ignore_index=True)  
                
        if 'lasair' in self.Brokers.names:
            idx = self.Brokers.names.index('lasair')
            res = self.Brokers.queryresults[idx]
            if len(res) == 0:
                print("query result is empty")
            else:
                if 'use_sherlock' in kwargs.keys() and not kwargs['use_sherlock']:
                    print("Lasair: no sherlock classification/crossmatch queried, use sdss/ned/glade for crossmatch")
                    idx = self.Brokers.names.index('lasair')
                    sncoor = self.Brokers.queryresults[idx][['oid','meanra','meandec']]
                    #crossmatch with sdss/glade first
                    dblist = [x for x in ['sdss','glade'] if x in self.CrossMatch.db_dict.keys()]
                    self.CrossMatch.run(sncoor=sncoor,select_db=dblist)        
                    xmatch_res = self.CrossMatch.result.sort_values('distance')

                    if 'ned' in self.CrossMatch.db_dict.keys():
                        #crossmatch ned for photoz in range
                        search_idx = (xmatch_res['photoz'] > 0.01) & (xmatch_res['photoz'] < 0.08) & (xmatch_res['specz'].isna())
                        if np.sum(search_idx) > 0:
                            print("Cross match NED for 0.01 < sdss_photoz < 0.08, num = ",np.sum(search_idx))
                            search_oid = xmatch_res.loc[search_idx,'oid'].drop_duplicates()
                            coor = sncoor.set_index('oid').loc[search_oid].reset_index()
                            self.CrossMatch.run(sncoor=coor,select_db=['ned'])
                            xmatch_res = pd.concat([xmatch_res,self.CrossMatch.result.sort_values('distance')],sort=False,ignore_index=True)

                        #crossmatch ned for dec<0
                        south_idx = sncoor['meandec'] < 0
                        if np.sum(south_idx) > 0:
                            print("Cross match NED for dec < 0, num = ",np.sum(south_idx))
                            coor = sncoor.loc[south_idx].drop_duplicates('oid')
                            self.CrossMatch.run(sncoor=coor,select_db=['ned'])
                            xmatch_res = pd.concat([xmatch_res,self.CrossMatch.result.sort_values('distance')],sort=False,ignore_index=True)

                    xmatch_res = xmatch_res.sort_values(['oid','distance'])   
                    res = self.Brokers.queryresults[idx].merge(xmatch_res,how='inner',on='oid')
                
                    res = self._assign_redshift_from_sdss(res)
                else:
                    res = self._assign_redshift_lasair(res)
                    
                res = self._rename_columns(res,orignames=['latestgmag','latestrmag'],
                                           newnames=['gmaglatest','rmaglatest'])
                res['Broker'] = 'lasair'
                self.results = pd.concat([self.results,res],sort=False,join='outer',ignore_index=True)
                
        if not realtime:
            from broker_query import Alerce
            alerce = Alerce()
            alerce.setup()
    
            print("In non-realtime mode querying detections for latest mag in mjd range.")
            
            oidlist = self.results.oid.unique()
            # query detection
            query="select oid, fid, mjd, magpsf, sigmapsf, rb from detections where oid in {} and mjd < {}".format(tuple(oidlist),mjdend)
            alerce.load_query(query)
            alerce.make_query()
            qresults = alerce.queryresult
            res = self._get_latestmag_before_mjd(qresults)
#             maglims = ((res.gmaglatest > 16) & (res.gmaglatest < 20)) | ((res.rmaglatest > 16) & (res.rmaglatest < 20))
            param_names = ['gmax','rmax','gmin','rmin']
            pdict = {}
            for p in param_names:
                if p in kwargs.keys():
                    pdict[p] = kwargs[p]
                else:
                    not_kwargs = [x for x in param_names if x not in kwargs.keys()]
                    pdict.update(self._get_default_parameter_value(self.Brokers.brokers[0].gen_querystr,param_names=not_kwargs))
#             print(pdict)
            maglims = ((res.gmaglatest > pdict['gmin']) & (res.gmaglatest < pdict['gmax'])) | \
                      ((res.rmaglatest > pdict['gmin']) & (res.rmaglatest < pdict['gmax']))
            oidres = res[maglims].oid.unique()
#             self.results = self.results[self.results.oid.isin(oidres)]
            cols = [x for x in self.results.columns if x not in ['gmaglatest','rmaglatest']]
            self.results = self.results[cols].merge(res[maglims][['oid','gmaglatest','rmaglatest']],on='oid',how='inner')
    
        if len(self.results) > 0:
            self.results.fillna(value=pd.np.nan, inplace=True)
            self.MakeCuts.load_data(self.results)            
            self.MakeCuts.makecuts(cuts=self.selection_cuts,**kwargs)
#             for flag in self.selection_cuts:
#                 cut_to_run = getattr(self.MakeCuts,flag)
#                 cut_to_run(**kwargs)
#                 self.MakeCuts.show_stats(flag='flag_'+flag)
            self.MakeCuts.plot(**kwargs)
            self.results = self.MakeCuts.data
        else:
            raise RuntimeError("Empty queries.")
                
    def _get_latestmag_before_mjd(self,data):
        func = lambda x: self.__get_latestmag_before_mjd_single_oid(x)
        res = data.groupby('oid').apply(func)
        data = data.merge(res,on='oid')
        return data
        
    def __get_latestmag_before_mjd_single_oid(self,detections):
        if len(detections['oid'].unique()) > 1:
            raise RuntimeError('Only single oid is allowed.')
        gmags = detections[detections['fid'] == 1]
        rmags = detections[detections['fid'] == 2]
        if len(gmags) > 0:
            gmaglatest = gmags.loc[gmags.mjd.idxmax()].magpsf
        else:
            gmaglatest = np.nan
        if len(rmags) > 0:
            rmaglatest = rmags.loc[rmags.mjd.idxmax()].magpsf
        else:
            rmaglatest = np.nan
        return pd.DataFrame({'gmaglatest':gmaglatest,'rmaglatest':rmaglatest},index=[0])
                
    def _assign_redshift_from_sdss(self,res):
        res['z'] = -99.
        res['zerr'] = -99.
        res.loc[res.specz==False,'specz'] = -99.
        res.loc[~res['specz'].isna(),'z'] = res.loc[~res['specz'].isna(),'specz']
        res.loc[~res['specz'].isna(),'zerr'] = 0.001
        idx = (~res['photoz'].isna()) & (res['z'] == -99.)
        res.loc[idx,'z'] = res.loc[idx,'photoz']
        res.loc[idx,'zerr'] = res.loc[idx,'photoz_err']
        return res
    
    def _assign_redshift_lasair(self,res):
        res['zerr'] = 0.001
        idx = res['z'].isna()
        res.loc[idx,'z'] = res.loc[idx,'photoz']
        res.loc[idx,'zerr'] = res.loc[idx,'photoz_err']
        return res
    
    def _rename_columns(self,res,orignames=[],newnames=[]):
        if len(orignames) != len(newnames):
            raise ValueError("newnames must have the same length as orignames")
        if len(orignames) > 0:
            for name1,name2 in zip(orignames,newnames):
                res[name2] = res[name1]
        return res
    
    def _get_default_parameter_value(self,function,param_names=None):
        res_dict = {}
        for p in param_names:
            v = inspect.signature(function).parameters[p].default
            res_dict[p] = v
        return res_dict    
        
class PipePro():
    def __init__(self):
        pass
    
    def run(self):
        pass
    
    
class Brokers(PipePro):       
    def __init__(self,brokers=[]):
        from broker_query import Alerce, Antares, Lasair, TNS
        self.names = brokers
        self.brokers = []
        for name in brokers:
            if name.startswith('alerce'):
                broker = Alerce()
            elif name.startswith('antares'):
                broker = Antares()
            elif name.startswith('lasair'):
                broker = Lasair()
            elif name.startswith('tns'):
                broker = TNS()
            self.brokers.append(broker)

    def run(self,realtime=True,mjdstart=None,mjdend=None,**kwargs):
        self.queryresults = []
        for broker in self.brokers:
            broker.setup()
            broker.gen_querystr(realtime=realtime,mjdstart=mjdstart,mjdend=mjdend,**kwargs)
            broker.make_query()
            self.queryresults.append(broker.queryresult)
    
class CrossMatch(PipePro):
    def __init__(self,catalogues=[]):
        self.db_dict = {}
        if len(catalogues) > 0:
            if 'sdss' in catalogues:
                from cross_match_host import SDSSQuery
                self.db_dict['sdss'] = SDSSQuery()
            if 'ned' in catalogues:
                from cross_match_host import NEDQuery
                self.db_dict['ned'] = NEDQuery()
            if 'glade' in catalogues:
                from cross_match_host import GLADEQuery
                self.db_dict['glade'] = GLADEQuery()
    
    def run(self,sncoor=None,select_db=['sdss']):
        result_df = pd.DataFrame()
        for key in self.db_dict.keys():
            if key in select_db:
                db = self.db_dict[key]
                db.run(sncoor)
                db.queryresult['xmatch_db'] = key
                result_df = pd.concat([result_df,db.queryresult],sort=False,ignore_index=True)
#         print(result_df)
        self.result = result_df.sort_values('oid')

class MakeCuts(PipePro):
    def __init__(self):
        pass

    def load_data(self,data):
        self.data = data
    
    def zcut(self,zlow=0.016,zhigh=0.08, zerrlim=0.02, **kwargs):
        print("Selecting {:.3f} < z < {:.3f}".format(zlow,zhigh))
        res = self.data
        res['flag_zcut'] = False
        idx = (res['z']>zlow) & (res['z']<zhigh)
        if zerrlim is not None:
            print(' and zerr < {:.3f}'.format(zerrlim))
            idx = idx & (res['zerr'] < zerrlim) & (res['zerr'] > 0.)
        res.loc[idx,'flag_zcut'] = True
        self.data = res
        
    def magzcut(self,dmag_min=1.,dmag_max=4.5,magabs=-19.1,**kwargs):
        from utils import cal_maglim
        res = self.data
        print("Selecting candidates based on mag vs z: {:.2f} < dmag (from max) < {:.2f}".format(dmag_min,dmag_max))
        res['flag_magzcut'] = False
        func = lambda x: cal_maglim(x,magabs=magabs)
        res['maglim'] = res['z'].apply(func)
        res['dmag_g'] = res['gmaglatest'] - res['maglim']
        res['dmag_r'] = res['rmaglatest'] - res['maglim']
        idx = ((res['dmag_g'] > dmag_min) & (res['dmag_g'] < dmag_max)) | ((res['dmag_r'] > dmag_min) & (res['dmag_r'] < dmag_max))
        res.loc[idx,'flag_magzcut'] = True
        self.data = res
        
    def rbcut(self,rb=0.5,**kwargs):
        
        from broker_query import Alerce
        alerce = Alerce()
        alerce.setup()
    
        print("Selecting candidates w/ rb >{:.2f}".format(rb))

        res = self.data
        if 'flag_rbcut' in res:
            res = res.drop(columns=['flag_rbcut'])
        oidlist = self.data.oid.unique()
        # query rb
        query="select oid, candid, rb from detections where oid in {}".format(tuple(oidlist))
        alerce.load_query(query)
        alerce.make_query()
        qresults = alerce.queryresult
        if len(qresults)> 0:
            rbcount = qresults.groupby('oid').apply(lambda x: (x.rb > rb).sum()) > 0
            rbcount.name = 'flag_rbcut'
            res = res.merge(rbcount,on='oid',how='outer')
            res.loc[res['flag_rbcut'].isna(),'flag_rbcut'] = True
        else:
            res['flag_rbcut'] = True
        res['flag_rbcut'] = res['flag_rbcut'].astype(bool)
        self.data = res
        
    def hostdistcut(self,dist_min=2,dist_max=None,**kwargs):
        """
        cut distance to host (arcsec)
        """
        print("Cut on distance to host: {} < dist < {} (Arcsec)".format(dist_min,dist_max))
        res = self.data
        res['flag_hostdistcut'] = False
        if 'distance' in res.columns and 'separationArcsec' in res.columns:
            idx = (res.distance > dist_min/60.) | (res.separationArcsec > dist_min)
            if dist_max is not None:
                idx = idx & ((res.distance < dist_max/60.) | (res.separationArcsec < dist_max))
        elif 'distance' in res.columns and 'separationArcsec' not in res.columns:
            idx = (res.distance > dist_min/60.)
            if dist_max is not None:
                idx = idx & (res.distance < dist_max/60.) 
        elif 'distance' not in res.columns and 'separationArcsec' in res.columns:
            idx = (res.separationArcsec > dist_min)
            if dist_max is not None:
                idx = idx & (res.separationArcsec < dist_max)
        res.loc[idx,'flag_hostdistcut'] = True
        self.data = res
        
    def ztflimcut(self,maglim=19.8,nobs=1,**kwargs):
        print("Cut on magnitude lim for nobs <= {}: maglim = {}".format(nobs,maglim))
        res = self.data
        res['flag_ztflimcut'] = False
#         idx = (res['nobs'] > nobs) | ((res['gmaglatest'] < maglim) | (res['gmaglatest'].isna())) & ((res['rmaglatest'] < maglim) | (res['rmaglatest'].isna()))
        idx = (res['nobs'] > nobs) | (res['gmaglatest'] < maglim) | (res['rmaglatest'] < maglim) | (res['Broker'] == 'tns')
        res.loc[idx,'flag_ztflimcut'] = True
        self.data = res
        
    def olddetectioncut(self,dayslim=5,fromday=None,**kwargs):
        from astropy.time import Time
        from datetime import datetime,date
        
        if fromday is None:
            fromday = Time(datetime.today()).mjd 

        print("Cut on days since last detection from {}: dt < {}".format(Time(fromday,format='mjd').iso,dayslim))
        res = self.data
        res['flag_olddetectioncut'] = False
        idx = fromday - res['lastmjd'] > dayslim
        res.loc[~idx,'flag_olddetectioncut'] = True
        self.data = res
        
    def visibility(self,startdate=None,enddate=None):
        pass
    
    def makecuts(self,cuts=[],show_stats=True,**kwargs):
        for flag in cuts:
            cut_to_run = getattr(self,flag)
            cut_to_run(**kwargs)
            if show_stats:
                self.show_stats(flag='flag_'+flag)
    
    def show_stats(self,flag=None):
        if flag is not None:
            ncut = np.sum(~self.data[flag])
            print("Number fails cut [{}]: {}/{}".format(flag,ncut,len(self.data)))
    
    def aftercuts(self):
        res = self.data
        flags = [x for x in res.columns if x.startswith('flag')]
        idx = np.array([True] * len(res))
        for flag in flags:
            idx = idx & res[flag]
        return res[idx]
    
    def plot(self,magz_plot=True,magz_band='r',magabs=-19.1,**kwargs):
        from utils import plot_maglims
        print("Plotting MakeCuts...")
#         plotdata = self.data[self.data.flag_zcut & self.data.flag_magzcut & self.data.flag_rbcut]
        plotdata = self.aftercuts()
        if 'zerr' in self.data.columns:
            plt.errorbar(plotdata['z'],plotdata['gmaglatest'],xerr=plotdata['zerr'],color='g',fmt='o',label='g')
            plt.errorbar(plotdata['z'],plotdata['rmaglatest'],xerr=plotdata['zerr'],color='r',fmt='o',label='r')
        else:
            plt.plot(plotdata['z'],plotdata['gmaglatest'],xerr=plotdata['zerr'],color='g',marker='o',label='g')
            plt.plot(plotdata['z'],plotdata['rmaglatest'],xerr=plotdata['zerr'],color='r',marker='o',label='r')
        plt.ylim(plt.ylim()[::-1])
        plt.legend()
        if magz_plot:
            plot_maglims(band=magz_band,magabs=magabs,**kwargs)
            
class EarlyClass(PipePro):
    def __init__(self):
        pass

