import re
import requests
from bs4 import BeautifulSoup
import json
from collections import OrderedDict
from io import StringIO
import pandas as pd
from astropy.time import Time
from datetime import datetime,date,timedelta
from tns_api_search import search, get, format_to_json, get_file
from astropy.coordinates import SkyCoord
from astropy import units as u

url_tns_api="https://wis-tns.weizmann.ac.il/api/get"
from credentials import tns
api_key = tns.settings.API_KEY

def goodrow(class_):
    return ((class_=="row-even public odd") or (class_=='row-odd public odd'))

def getTNS(reportdays=5,discoverdays=5,enddate=None,classified=1,disc_mag_min=16,disc_mag_max=21,z_min=0.015,z_max=0.08,
           skip_ztf=False,num_page=100,verbose=False,otherparams={},**kwargs):
    '''
    returns a coma separated list with the redshift, internal name, discovery date, and discovery magnitude
    of objetcs from TNS which match the search criteria

    parameters:
        reportdays - maximum number of days that have past since being reported
        z_min - minimum redshift
        z_max - maximum redshift
        disc_mag_min - minimum discovery magnitude
        disc_mag_max - maximum discovery magnitude
            Note: I believe this is just a numerical cut, not physical, i.e. the minimum is the lowest numerical
                    value that will be returned
        calssified - 1: is classified, 0: classification not considered
        unclassified - 1: is unclassified, 0: unclassification not considered
    '''
#     link = f'https://wis-tns.weizmann.ac.il/search?&discovered_period_value={reportdays}&discovered_period_units=days&unclassified_at={unclassified}&classified_sne={classified}&name=&name_like=0&isTNS_AT=all&public=all&coords_unit=arcsec&redshift_min={z_min}&redshift_max={z_max}&discovery_mag_min={disc_mag_min}&discovery_mag_max={disc_mag_max}&objtype=3&sort=desc&order=discoverydate&num_page=500'
    link = 'https://wis-tns.weizmann.ac.il/search'
    if enddate is None:
        enddate = date.today()
    startdate = enddate - timedelta(discoverdays)
    
    params = {"discovered_period_value":reportdays,
              "discovered_period_units":"days",
              "date_start[date]":startdate.isoformat(),
              "date_end[date]":enddate.isoformat(),
              "classified_sne":int(classified),
              "unclassified_at":int(not(classified)),
              "discovery_mag_min":disc_mag_min,
              "discovery_mag_max":disc_mag_max,
              "num_page":num_page
             }
    params.update(otherparams)
    if classified:
        params.update({"objtype":3,
                       "redshift_min":z_min,
                       "redshift_max":z_max,
                       "sort":"desc",
                       "order":"discoverydate"})
    else:
        params.update({"at_type":1,
                       "sort":"asc",
                       "order":"internal_name"})

    r = requests.get(link,params=params)
    if verbose:
        print(r.url)
    soup = BeautifulSoup(r.text, "lxml")
    
    return_arr = []
    tr = soup.find_all('tbody')
    if verbose:
        print("Number of tables on the webpage:",len(tr))
    if len(tr)>0:
        tr = tr[0]
    else:
        raise RuntimeError("No result is found")
               
    cols = ['internal_name','redshift','ra','decl','hostname','host_redshift','discoverydate','discoverymag','disc_filter_name','name','ot_name']
    dflist = []
    if verbose:
        print("Number of rows in search result: ",len(tr.find_all(class_=goodrow,recursive=False)))
    for row in tr.find_all(class_=goodrow,recursive=False):
        df = {}
        for col in cols:
            value = row.find('td',class_='cell-{}'.format(col),recursive=False)
            if value is None:
                df[col] = None
            else:
                df[col] = value.text
        df['name'] = df['name'].split()[1]
        if (not classified) & skip_ztf & df['internal_name'].startswith('ZTF'):
            break
        dflist.append(df)
            
    df = pd.DataFrame(dflist)
    df.columns = ['internal_name','redshift','ra_s','dec_s','hostname','host_redshift','discoverydate','discoverymag','disc_filter_name','tns_name','type']
    c = SkyCoord(ra=df.ra_s.values, dec=df.dec_s.values, unit=(u.hourangle,u.deg))
    df['meanra'] = c.ra.degree
    df['meandec'] = c.dec.degree
    df['oid'] = df['tns_name']
    return df.sort_values('discoverydate',ascending=False).reset_index(drop=True)


def get_tns_name(internal_name):
    search_obj = [("internal_name",internal_name)]
    response = search(url_tns_api,search_obj)
    if None not in response:
        json_data = format_to_json(response.text)
        reply = json.loads(json_data)['data']['reply']
        if len(reply) == 0:
            return
        else:
            return reply[0]['objname']
    else:
        print(response[1])
        return
        
def get_tns_data(tns_name):
    data = {}
    get_obj=[("objname",tns_name), ("photometry","1"), ("spectra","1")]
    response=get(url_tns_api,get_obj)
    if None not in response:
        # Here we just display the full json data as the response
        json_data = format_to_json(response.text)
        data['meta'] = format_meta(json.loads(json_data)['data']['reply'])
        photometry = json.loads(json_data)['data']['reply']['photometry']
        spectra = json.loads(json_data)['data']['reply']['spectra']
        data['photometry'] = format_photometry(photometry)
        data['spectra'] = format_spectra(spectra)
    else:
        print (response[1])
        data = None
    return data
        
def format_meta(reply):
    cols = ['internal_name','redshift','radeg','decdeg','hostname','host_redshift','discoverydate','discoverymag','discmagfilter','name']
    df = {k: reply[k] for k in cols}
    return pd.DataFrame([pd.DataFrame(df).loc['name']])
    
def format_photometry(photometry):
    dflist = []
    for epoch in photometry:
        dflist.append(epoch)
    df = pd.DataFrame(dflist)
    cols = ['flux_unit','instrument','telescope','filters']
    for col in cols:
        df[col] = df[col].apply(pd.Series)['name']
    df['mjd'] = Time(list(df['obsdate'].values),format='iso').mjd       
    return df

def format_spectra(spectra):
    if len(spectra) == 0:
        return
    else:
        dflist = []
        for epoch in spectra:
            dflist.append(epoch)
        df = pd.DataFrame(dflist)
        cols = ['source_group', 'instrument', 'telescope']
        for col in cols:
            if col in df.columns:
                df[col] = df[col].apply(pd.Series)['name']
        return df
    