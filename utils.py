import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
# os.environ['SFD_DIR'] = '/home/mi/sfddata-master'
# os.environ['SIRAHPIPE_DIR'] = '/home/mi/ztf/sirah_target_pipe'
import sncosmo
from scipy.interpolate import interp1d
from astropy.time import Time
from datetime import datetime,date
from io import BytesIO
import pickle
from matplotlib.backends.backend_pgf import PdfPages
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy.coordinates import SkyCoord
from tnsScrape import get_tns_name
from astropy.table import Table

def get_mu(z):
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    return cosmo.distmod(z).value

def cal_maglim(z, magabs=-19.1):
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    z = np.array(z)
    mag = magabs + cosmo.distmod(z).value
    return mag
    
def plot_maglims(band='r',z=0.,phase=[0.,-4.,-7.,-10.,-13.,-15.],magabs=-19.1):
    zarr = np.linspace(0.01,0.08,30)
    mag = cal_maglim(zarr, magabs=magabs)
    dmag = get_dmag(band,z,phase=phase)['dmag']
    for d,p in zip(dmag,phase):
        plt.plot(zarr, mag+d, label='phase={:.1f} (dmag={:.2f})'.format(p,d))  
    plt.axhline(y=20.5,ls=':',color='k',label='ztf limit')
    plt.legend(bbox_to_anchor=(1.05,1))
    plt.xlabel('z')
    plt.ylabel('mag')


def query_timeseries(oid,plotlc=False,source='ztf',broker='alerce'):
    if source == 'ztf':
        return _query_timeseries_ztf(oid,plotlc=plotlc,broker=broker)
    elif source == 'tns':
        return _query_timeseries_tns(oid,plotlc=plotlc)
    else:
        return
    
def _query_timeseries_tns(oid,plotlc=False):
    
    from tnsScrape import get_tns_data
    data = get_tns_data(oid)
    phot = data['photometry']
    meta = data['meta']
    ra = meta['radeg']
    dec = meta['decdeg']
    res = {}
    idx_nondet = phot.flux.isna()
    res['LC_nondet'] = phot[idx_nondet]
    res['LC_det'] = phot[~idx_nondet]
    func = lambda row: mwebv_corr(ra,dec,'g')
    res['LC_det']['mag_dered'] = res["LC_det"]['flux'] - res["LC_det"].apply(func,axis=1)    
    idx = np.array([x in 'gr' for x in res['LC_det']['filters']])
    if np.sum(idx)>0:
        func = lambda row: mwebv_corr(ra,dec,row['filters'])
        res["LC_det"].loc[idx,'mag_dered'] = res["LC_det"].loc[idx,'flux'] - res["LC_det"].loc[idx].apply(func,axis=1)    

    return res
    
def _query_timeseries_ztf(oid,plotlc=False,broker='alerce'):
    """
    adopted from alerce notebooks
    """
    if broker == 'alerce':
        from broker_query import Alerce
        alerce = Alerce()
        alerce.setup()

        results = {}

        # query detections and sort by mjd
        query="select oid, candid, ra, dec, fid, mjd, magpsf, sigmapsf,rb from detections where oid='%s'" % oid
        alerce.load_query(query)
        alerce.make_query(verbose=False)
        results["LC_det"] = alerce.queryresult
        labels = {1: 'g', 2: 'r'}
        func = lambda row: mwebv_corr(row['ra'],row['dec'],labels[row['fid']])
        results["LC_det"]['magpsf_dered'] = results["LC_det"]['magpsf'] - results["LC_det"].apply(func,axis=1)    

        # query non detections and sort by mjd
        query="select oid, fid, mjd, diffmaglim from non_detections where oid='%s'" % oid
        alerce.load_query(query)
        alerce.make_query(verbose=False)
        results["LC_nondet"] = alerce.queryresult
        
    elif broker == 'lasair':
        from broker_query import Lasair
        lasair = Lasair()
        lasair.setup()
        results = {}

        # query detections and sort by mjd
        query='''SELECT candidates.objectId as oid, candidates.candid, candidates.ra, candidates.decl, candidates.fid, 
                 candidates.jd-2400000.5 as mjd, candidates.magpsf, candidates.sigmapsf, candidates.drb as rb 
                 FROM candidates WHERE candidates.objectId = "%s"''' % oid
        lasair.load_query(query)
        lasair.qcols = ['oid','candid','ra','dec','fid','mjd','magpsf','sigmapsf','rb']
        lasair.make_query(verbose=False)
        results["LC_det"] = lasair.queryresult
        labels = {1: 'g', 2: 'r'}
        func = lambda row: mwebv_corr(row['ra'],row['dec'],labels[row['fid']])
        results["LC_det"]['magpsf_dered'] = results["LC_det"]['magpsf'] - results["LC_det"].apply(func,axis=1)    

        # query non detections and sort by mjd
        query='''select noncandidates.objectId, noncandidates.fid, noncandidates.jd-2400000.5 as mjd, noncandidates.diffmaglim 
               from noncandidates where noncandidates.objectId="%s"''' % oid
        lasair.load_query(query)
        lasair.qcols = ['oid','fid','mjd','diffmaglim']
        lasair.make_query(verbose=False)
        results["LC_nondet"] = lasair.queryresult

    if plotlc:
        plotLC(oid, results["LC_det"], results["LC_nondet"], source='ztf')
        plt.show()
        
    return results

def get_dmag(band,z,phase=None,x1=0.):
    register_ztf_bandpass(band)
    if phase is None:
        phase = np.arange(-17,20,1)
    else:
        phase = np.array(phase)
    model = sncosmo.Model(source='salt2')
    model.set(z=z)
    model.set(x1=x1)
    time = phase*(1.+z)
    modelmag = model.bandmag('ztf'+band,'ab',time) #bandmag() takes time instead of phase
    modelmag_peak = model.bandmag('ztf'+band,'ab',0.)
    dmag = modelmag - modelmag_peak  
    return {'dmag':dmag,'time':time,'z':z} 

def interp_lc(t,lc,mag=0.):
    idx = np.isfinite(lc)
    interpf = interp1d(t[idx],lc[idx],fill_value='extrapolate') 
    func = lambda x: interpf(x) - mag   
    return func

def get_peakmjd(func,mjd):
    from scipy.optimize import root_scalar
    rootres = root_scalar(func,x0=0.,method='bisect',bracket=[-20,0])
    p = rootres.root      
    peakmjd = mjd - p     
    return peakmjd

def plot_lc_prediction(band,z,mag,mjd,magabs=-19.1):
    
    dmagres = get_dmag(band,z, phase=np.linspace(-17,20,50))
    dmag = dmagres['dmag']
    dmagres_pos_x1 = get_dmag(band,z, phase=np.linspace(-17,20,50),x1=2.)
    dmag_pos = dmagres_pos_x1['dmag']
    dmagres_neg_x1 = get_dmag(band,z, phase=np.linspace(-17,20,50),x1=-2.)
    dmag_neg = dmagres_neg_x1['dmag']

    t = dmagres['time']
    peakmag = cal_maglim(z,magabs=magabs)
    lc = dmag + peakmag 
    lc_pos = dmag_pos + peakmag - 0.2
    lc_neg = dmag_neg + peakmag + 0.2
    func = interp_lc(t,lc,mag=mag)
    func_pos = interp_lc(t,lc_pos,mag=mag)
    func_neg = interp_lc(t,lc_neg,mag=mag)
    
    result = {}

    try:
        peakmjd = get_peakmjd(func,mjd)
        time = t + peakmjd
        lc_plot = interp_lc(time,lc)(time)
        plt.plot(time,lc_plot,color=band,ls='--')
        plt.axvline(x = peakmjd, label='estimated peak in {}'.format(band),color=band,ls=':')
        result['peakmjd'] = peakmjd
    except:
        print("Can't estimate phase")
        if func(-20)*func(0)>0:
            print("The given mag is out of the salt2 model range given the redshift")
            result['peakmjd'] = np.nan
    try:
        peakmjd_pos = get_peakmjd(func_pos,mjd)
        time_pos = t + peakmjd_pos
        lc_pos_plot = interp_lc(time_pos,lc_pos)(time)
        plt.fill_between(time, lc, lc_pos_plot, color=band, alpha=0.2)
        plt.axvspan(peakmjd, peakmjd_pos, alpha=0.2, color=band)
    except:
        print("Can't estimate phase")
        if func_pos(-20)*func_pos(0)>0:
            print("Higer Limit ({}): The given mag is out of the salt2 model range given the redshift".format(band))
    try:        
        peakmjd_neg = get_peakmjd(func_neg,mjd)
        time_neg = t + peakmjd_neg
        lc_neg_plot = interp_lc(time_neg,lc_neg)(time)
        plt.fill_between(time, lc, lc_neg_plot, color=band, alpha=0.2)
        plt.axvspan(peakmjd, peakmjd_neg, alpha=0.2, color=band)
    except:
        print("Can't estimate phase")
        if func_neg(-20)*func_neg(0)>0:
            print("Lower Limit ({}): The given mag is out of the salt2 model range given the redshift".format(band))
    
    return result

def get_ps1_image(ra,dec,size=240):
    """
    size: pixel number for 0.25 arcsec/pixel
    """
    from PIL import Image
    import requests
    from ps1_image import geturl

    try:
        url = geturl(ra, dec, size=size, output_size=None, filters="grz", format="jpg", color=True)
        r = requests.get(url)
        im = Image.open(BytesIO(r.content))
    except:
        print("Can't get ps1 image")
        im = None
    return im
    
def plot_aladin(ra,dec):
    from ipywidgets import Layout
    aladinres = aladin(ra,dec)
    aladinres.layout = Layout(width='40%')
    display(aladinres)

def aladin(ra,dec):
    import ipyaladin as ipyal
    target='{} {}'.format(ra,dec)
    aladinres = ipyal.Aladin(target=target,fov=0.02,survey='P/PanSTARRS/DR1/color/z/zg/g')
    return aladinres

def plotLC(oid, SN_det, SN_nondet, source='ztf', figsize=None, plot_ylim=(21,15)):
    if source == 'ztf':
        return _plotLC_ztf(oid, SN_det, SN_nondet, figsize=figsize, plot_ylim=plot_ylim)
    elif source == 'tns':
        return _plotLC_tns(oid, SN_det, SN_nondet, figsize=figsize, plot_ylim=plot_ylim)
    else:
        return

def _plotLC_tns(oid, SN_det, SN_nondet, figsize=None,plot_ylim=(21,15)):

    fig, ax = plt.subplots(figsize = figsize)
    
    SN_det = SN_det.sort_values('mjd')
    if len(SN_nondet) > 0:
        SN_nondet = SN_nondet.sort_values('mjd')   
    SN_det.loc[SN_det.fluxerr.isna(),'fluxerr'] = 0.
    
    colorlist = 'bmyk'
    i=0
    for f in SN_det.filters:
        if f in ['cyan','orange','g','r']:
            color = f
        else:
            color = colorlist[i]
            i += 1
        mask = SN_det.filters == f
        if np.sum(mask) > 0:
            ax.errorbar(SN_det[mask].mjd, SN_det[mask].flux, yerr=SN_det[mask].fluxerr,
                        marker = '*', ls=':', label = f, color=color,alpha=0.6)
            ax.errorbar(SN_det[mask].mjd, SN_det[mask].mag_dered, yerr=SN_det[mask].fluxerr,
                        marker = 'o', label = f+'_corr', color=color)

    if len(SN_nondet) > 0:
        for f in SN_nondet.filters:
            if f in ['cyan','orange','g','r']:
                color = f
            else:
                color = None
            mask = SN_nondet.filters == f
            if np.sum(mask) > 0:
                ax.scatter(SN_nondet[mask].mjd, SN_nondet[mask].limflux,
                           marker = 'v', ls=':', label = "lim.mag %s"%f, color=color,alpha=0.5) 
            
    ax.set_title(oid)
    ax.set_xlabel("MJD")
    ax.set_ylabel("Magnitude")
    ax.legend()
    ax.set_ylim(plot_ylim)
    ax.set_xlim((SN_det.mjd.max()-30,SN_det.mjd.max()+30))
    
    return fig, ax


def _plotLC_ztf(oid, SN_det, SN_nondet, figsize=None, plot_ylim=(21,15)):
    """
    adopted from alerce notebooks
    """
    fig, ax = plt.subplots(figsize = figsize)
    labels = {1: 'g', 2: 'r'}
    colors = {1: 'g', 2: 'r'}
    SN_det = SN_det.sort_values('mjd')
    if len(SN_nondet) > 0:
        SN_nondet = SN_nondet.sort_values('mjd')
    for fid in [1, 2]:
        mask = SN_det.fid == fid
        if np.sum(mask) > 0:            
            ax.errorbar(SN_det[mask].mjd, SN_det[mask].magpsf, 
                        yerr = SN_det[mask].sigmapsf, c = colors[fid], marker = '*', ls=':', label = labels[fid], alpha=0.6)
            ax.errorbar(SN_det[mask].mjd, SN_det[mask].magpsf_dered, 
                        yerr = SN_det[mask].sigmapsf, c = colors[fid], marker = 'o', label = labels[fid]+'_corr')
        if len(SN_nondet) > 0:
            mask = (SN_nondet.fid == fid) & (SN_nondet.diffmaglim > -900)
            if np.sum(mask) > 0:            
                ax.scatter(SN_nondet[mask].mjd, SN_nondet[mask].diffmaglim, c = colors[fid], alpha = 0.5,
                    marker = 'v', label = "lim.mag. %s" % labels[fid])
    ax.set_title(oid)
    ax.set_xlabel("MJD")
    ax.set_ylabel("Magnitude")
    ax.legend()
    ax.set_ylim(plot_ylim)
    ax.set_xlim((SN_det.mjd.max()-30,SN_det.mjd.max()+30))

    return fig, ax

def plot_extra_lc(SN, plot_prediction=False, z=None,magabs=None):
    colormap = {'o':'orange','c':'cyan','g':'g','r':'r','orange':'orange','cyan':'cyan'}
    for f in SN['filter'].unique():
        SN_f = SN[SN['filter']==f].sort_values('mjd')
        plt.errorbar(SN_f.mjd.values,SN_f.mag.values,yerr=SN_f.magerr.values,color=colormap[f], marker='o',mfc='white',label=f)
        if plot_prediction:
            mjd = SN[SN['filter'] == f].mjd.max()
            mag = SN[(SN['filter'] == f) & (SN.mjd == mjd)].mag.values[0]
            if f in 'gr':
                res = plot_lc_prediction(f,z,mag,mjd,magabs=magabs)
            elif f[0] in 'oc':
                fmap = {'o':'r','c':'g'}
                res = plot_lc_prediction(fmap[f[0]],z,mag,mjd,magabs=magabs)
            else:
                res = plot_lc_prediction('g',z,mag,mjd,magabs=magabs)
        else:
            res = {'peakmjd': np.inf}
    return res
    
def get_mwebv(ra,dec):
    import sfdmap
    
    dustmap = sfdmap.SFDMap(os.getenv('SFD_DIR'))
    mwebv = dustmap.ebv(ra,dec)
    return mwebv
    
def fit_salt2(SN_det,z,mwebv):
    model = sncosmo.Model(source='salt2',
                          effects=[sncosmo.F99Dust()],
                          effect_names=['mw'],
                          effect_frames=['obs'])
    model.set(z=z)
    model.set(mwebv=mwebv)
    filt_map = {1:'g',2:'r'}
    zp = 27.5
    SN_det['flux'] = np.power(10.,-0.4*(SN_det.magpsf-zp))
    SN_det['fluxerr'] = np.absolute(0.921*SN_det.flux*SN_det.sigmapsf)
    SN_det['zp'] = zp
    SN_det['zpsys'] = 'ab'
    SN_det['filter'] = ['ztf'+filt_map[x] for x in SN_det.fid]
    SN = Table.from_pandas(SN_det)
    res, fitmodel = sncosmo.fit_lc(SN,model,['t0','x0','x1','c'],
                                   bounds={'x0':(0,1.),
                                           'x1':(-5.,5.),
                                           'c':(-3.,3.)})
    return res

def plot_salt2(modelpars,band,phase=None):
    if not band in 'gr':
        return
    else:
        model = sncosmo.Model(source='salt2')
        model.update(modelpars)
        if phase is None:
            phase = np.linspace(-15,20,50)
            time = phase*(1. + modelpars['z']) + modelpars['t0']
            mag = model.bandmag('ztf'+band,'ab',time)
        plt.plot(time,mag,'-',color=band,lw=1.5)
    
def mwebv_corr(ra,dec,band):

    mwebv = get_mwebv(ra,dec)
    
    model_nomw = sncosmo.Model(source='salt2')
    model_mw = sncosmo.Model(source='salt2',
                             effects=[sncosmo.F99Dust()],
                             effect_names=['mw'],
                             effect_frames=['obs'])
    model_mw.set(mwebv=mwebv)
    try:
        bandpass = sncosmo.get_bandpass('ztf'+band)
    except:
        register_ztf_bandpass(band)
        bandpass = sncosmo.get_bandpass('ztf'+band)
    
    mag_nomw = model_nomw.bandmag(bandpass,'ab',0.)
    mag_mw = model_mw.bandmag(bandpass,'ab',0.)
    mwebv_corr = mag_mw - mag_nomw
    
    return mwebv_corr

def register_ztf_bandpass(band):
    bandtable = pd.read_csv(os.path.expandvars("$SIRAHPIPE_DIR/data/ztf_filters/ZTF_{}.csv".format(band)),names=['wave','trans'],comment='#')
    bandtable = bandtable.drop_duplicates('wave')
    bandpass = sncosmo.Bandpass(wave=bandtable.wave,trans=bandtable.trans,wave_unit=u.nm)
    bandpass.name = 'ztf'+band
    sncosmo.register(bandpass,force=True)
    
    
def gen_outside_links(oid,source=None):
    res = {}
    res['lasair'] = 'https://lasair.roe.ac.uk/object/{}'.format(oid)
    res['alerce'] = 'https://alerce.online/object/{}'.format(oid)
    res['tns'] = 'https://wis-tns.weizmann.ac.il/object/{}'.format(oid)
    res['yse'] = 'https://ziggy.ucolick.org/yse/transient_detail/{}'.format(oid)
    res['snex'] = 'https://supernova.exchange/view_object?name={}'.format(oid)
    if source in res.keys():
        return res[source]
    else:
        raise ValueError("source not found")
        
def gen_plots(row,interactive=False,pdf_file=None,magabs=-19.1,extra_lc=False,lc_prediction=True,update_lc_prediction=False,
              last_detection_max=5,source='ztf',plot_ylim=(21,15),broker='alerce',ps1_image_size=320,plot_salt2_lc=False):
    
    from IPython.core.display import HTML
    today = Time(datetime.today()).mjd 
    weekday_today = datetime.today().weekday()
    next_tuesday = today - weekday_today + 8 if weekday_today > 1 else today - weekday_today + 1
    
    oid = row['oid']
    z = row['z']
    zerr = row['zerr']
    print("{}: z={:.4f} +/- {:.4f}".format(oid, z, zerr))
    ra = row['meanra']
    dec = row['meandec']
    mwebv = get_mwebv(ra,dec)
    c = SkyCoord(ra, dec, unit='deg')
    ra_precision = 3
    ra_s = c.ra.to_string(u.hourangle,sep=':',pad=True,precision=ra_precision)
    dec_s = c.dec.to_string(u.degree,sep=':',alwayssign=True,pad=True,precision=ra_precision-1)
    print("ra = {} dec = {} mwebv={:.3f}".format(ra_s, dec_s, mwebv))
    
    figures = []
    
    links = []
    if row['Broker'] in ['lasair','alerce']:
        for s in ['lasair','alerce']:
            link = gen_outside_links(oid,source=s)
            links.append(link)
            display(HTML('<a href="{}">{}</a>'.format(link,link)))
        tns_name = get_tns_name(oid)
        if tns_name is not None:
            for s in ['tns','yse','snex']:
                link = gen_outside_links(tns_name,source=s)
                links.append(link)
                display(HTML('<a href="{}">{}</a>'.format(link,link)))
    elif row['Broker'] == 'tns':
        for s in ['tns','yse','snex']:
            link = gen_outside_links(oid,source=s)
            links.append(link)
            display(HTML('<a href="{}">{}</a>'.format(link,link)))
    res = query_timeseries(oid,plotlc=False,source=source,broker=broker)
#     print(res) 

    # light curve
    fig, ax = plotLC(oid, res["LC_det"], res["LC_nondet"],figsize=(7,5),source=source,plot_ylim=plot_ylim)
    lc_det = res['LC_det']
    
    if plot_salt2_lc:
        salt2res = fit_salt2(lc_det,z,mwebv)
        salt2pars = {}
        for pname,p in zip(salt2res['param_names'],salt2res['parameters']):
            salt2pars[pname] = p
        print(salt2pars)
        if 'mwebv' in salt2pars.keys():
            del salt2pars['mwebv']
    
    if source == 'ztf':
        if 1 in lc_det.fid.unique():
            mjd_g = lc_det[lc_det.fid == 1].mjd.max()
            mag_g = lc_det[(lc_det.fid == 1) & (lc_det.mjd == mjd_g)].magpsf_dered.values
            if lc_prediction:
                resg = plot_lc_prediction('g',z,mag_g,mjd_g,magabs=magabs)
            else:
                resg = {'peakmjd': np.inf}
            if plot_salt2_lc:
                plot_salt2(salt2pars,'g')
        else:
            mjd_g = -99.
            resg = {'peakmjd': np.inf}
        if 2 in lc_det.fid.unique():
            mjd_r = lc_det[lc_det.fid == 2].mjd.max()
            mag_r = lc_det[(lc_det.fid == 2) & (lc_det.mjd == mjd_r)].magpsf_dered.values
            if lc_prediction:
                resr = plot_lc_prediction('r',z,mag_r,mjd_r,magabs=magabs) 
            else:
                resr = {'peakmjd': np.inf}
            if plot_salt2_lc:
                plot_salt2(salt2pars,'r')
        else:
            mjd_r = -99.
            resr = {'peakmjd': np.inf}
            
    else:
        lc_det = lc_det.drop_duplicates()
        mjd_g = lc_det.mjd.max()
        mag_g = lc_det[lc_det.mjd == mjd_g].mag_dered.values
        if lc_prediction:
            resg = plot_lc_prediction('g',z,mag_g,mjd_g,magabs=magabs)
        else:
            resg = {'peakmjd': np.inf}
        mjd_r = -99.
        resr = {'peakmjd': np.inf}
    
    if extra_lc:
        lcfile = 'data/extra_photometry/{}.txt'.format(oid)
        SN = pd.read_csv(lcfile)
        if not update_lc_prediction:
            plot_extra_lc(SN)
        else:
            res_extra = plot_extra_lc(SN, plot_prediction=True, z=z, magabs=magabs)
            
    #exclude those whose latest dectections are old
    if today-mjd_g > last_detection_max and today-mjd_r > last_detection_max:
        plt.clf()
        print("Latest detection is > {} days".format(last_detection_max))
        return {'too_old': True}
        
    plt.axvline(x=today,color='k',ls=':',label='today')
    plt.axvline(x=next_tuesday,color='b',ls=':',label='next tuesday')
    plt.legend(bbox_to_anchor=(1.15, 1), loc='upper left', borderaxespad=0.)
    
    peakmjd = np.min([resg['peakmjd'],resr['peakmjd']])
    if update_lc_prediction and extra_lc:
        peakmjd = np.min([peakmjd,res_extra['peakmjd']])
    if plot_salt2_lc:
        peakmjd = salt2pars['t0']
        plt.axvline(x=peakmjd,color='purple',ls=':',label='salt2-t0')
        plt.legend(bbox_to_anchor=(1.15, 1), loc='upper left', borderaxespad=0.)
    
    ax2 = ax.twiny()
    newticks = np.arange(-50,20,5)
    ticks = newticks*(1.+z) + peakmjd
#     ax2.plot((lc_det.mjd - peakmjd)/(1.+z), np.zeros(len(lc_det.mjd)))
    ax2.set_xticks(ticks)
    ax2.set_xticklabels(newticks)
    ax2.set_xlim(ax.get_xlim())
    
    ax3 = ax.twinx()
    mu = get_mu(z)
    newticks = np.arange(-25,-10,1)
    ticks = newticks + mu
    ax3.set_yticks(ticks)
    ax3.set_yticklabels(['%d'%int(t) for t in newticks])
    ax3.set_ylim(ax.get_ylim())
    for i in range(len(links),5):
        links.append('')
    
    if pdf_file:
        title = """{}: z={:.4f} +/- {:.4f}\n ra={}, dec={}, mwebv={:.3f}\n\n""".format(oid, z, zerr, ra_s, dec_s,mwebv)
        title += '\n'.join(links)
        title += '\n'
        ax.set_title(title)
        plt.tight_layout()
        figures.append(fig)
    plt.show()
        
    # mag/z
    fig = plt.figure(figsize=(7,3))
    plt.subplot()
    if source == 'ztf':
        if 1 in lc_det.fid.unique():
            plt.errorbar(z,mag_g,xerr=zerr,fmt='g*',label='g-latest (offset from tuesday={:.1f} days)'.format(next_tuesday-mjd_g))
        if 2 in lc_det.fid.unique():
            plt.errorbar(z,mag_r,xerr=zerr,fmt='r*',label='r-latest (offset from tuesday={:.1f} days)'.format(next_tuesday-mjd_r))
    else:        
        plt.errorbar(z,mag_g,xerr=zerr,fmt='g*',label='mag-latest (offset from tuesday={:.1f} days)'.format(next_tuesday-mjd_g))

    plt.ylim(plot_ylim)
    plot_maglims(band='r',z=z,magabs=magabs)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.title(oid)
    plt.tight_layout()
    plt.show()
    figures.append(fig)
    
    # host image
    if interactive:
        plot_aladin(ra,dec)
    else:
        fig,ax = plt.subplots(figsize=(7,3))
        im = get_ps1_image(ra,dec,size=ps1_image_size)
        if im is not None:
            ax.imshow(im)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.scatter(np.average(plt.xlim()),np.average(plt.ylim()),marker='+',color='yellow')
            ax.set_title(oid)
            figures.append(fig)
            plt.show()
        else:
            print("Failed to get image")
        
    if pdf_file:
        with PdfPages(pdf_file) as pdf:
            for fig in figures:
                pdf.savefig(fig)
        
    return {'too_old': False, 'phase_tuesday': next_tuesday-peakmjd}