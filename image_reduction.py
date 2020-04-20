from photutils import SkyCircularAperture,aperture_photometry
from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import pandas as pd
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground
from astropy.io import fits
        
class ImageData():
    def __init__(self,objectname=None,image_file=None,se_file=None):
        self.objectname = objectname
        self.image_file = image_file
        self.se_file = se_file
    
    def photometry(self,target,references,method='sextracter',**kwargs):
        if method.startswith('se'):
            p = Sextracter(self.se_file)
        elif method.startswith('ap'):
            p = AperturePhotometry(self.image_file)
        res = p.run(target,references,**kwargs)
        self.result = res
    
    def get_image_info(self):
        hdul = fits.open(self.image_file)
        image_info = hdul[1].header
        self.info = image_info
    
    def to_file(self,fname,mwebv_corr=0.):
        self.get_image_info()
        mjd = self.info['MJD-OBS']
        flt = self.info['FILTER'][0]
        self.result.update({'filter':flt,'mjd':mjd,'name':self.objectname})
        self.result['mag'] = self.result['mag'] - mwebv_corr
        df = pd.DataFrame(self.result,index=[0])
        df.to_csv(fname)
    
class AperturePhotometry():
    def __init__(self,image_file=None):
        if image_file is not None:
            self.load_image(image_file)
            self.get_background()
    
    def load_image(self,image_file):
        self.image_file = image_file
        hdul = fits.open(self.image_file)
        self.data = hdul[1].data
        self.wcs = WCS(hdul[1].header) 
        
    def get_background(self):
        print("Calculating background...")
        self.background = Background2D(self.data, (1024, 1024))
        print("Done.")
        
    def make_aperture(self,ra,dec,unit=u.deg,aperture_radius=1.):
        position = SkyCoord(ra,dec, unit=unit, frame='icrs')
        ap = SkyCircularAperture(position, r=aperture_radius * u.arcsec)
        return ap
        
    def get_mag(self,ap):
        tbl = aperture_photometry(self.data-self.background.background,ap,wcs=self.wcs)
        print(tbl)
        result = tbl['aperture_sum'][0]
        mag = -2.5*np.log10(result)
        return mag
    
    def run(self,target,references,aperture_radius=1.):
        target_ap = self.make_aperture(target['ra'],target['dec'],aperture_radius=aperture_radius)
        target_mag_ap = self.get_mag(target_ap)

        maglist = []
        
        for i,ref in enumerate(references):
            ref_ap = self.make_aperture(ref['ra'],ref['dec'],aperture_radius=aperture_radius)
            ref_mag_ap = self.get_mag(ref_ap)
            magi = target_mag_ap - ref_mag_ap + ref['obsmag']
            maglist.append(magi)
        mag = np.average(maglist)
        magerr = np.std(maglist) if i>0 else None          

        return {'mag':mag,'magerr':magerr}    
    
class Sextracter():
    def __init__(self,se_file=None):
        if se_file is not None:
            self.datatable = self.read_se_result(se_file)
    
    def read_se_result(self,se_file):
        res = ascii.read(se_file,format='sextractor')
        return res
        
    def match_object(self,position,table,match_cols={'ra':'ALPHA_J2000','dec':'DELTA_J2000'},tolerance=1e-3):
        index = np.array([True]*len(table))
        for key in ['ra','dec']:
            index = index & np.isclose(table[match_cols[key]],position[key],atol=tolerance)
        if np.sum(index) == 0:
            raise RuntimeError("No matched object given the position")
        elif np.sum(index) > 1:
            print(table[index])
            if tolerance < 1e-5:
                raise RuntimeError("More than one matched objects are found at given position. Try decreasing tolerance.")
            else:
                res = table[index][0:1]
                return res
        else:
            res = table[index]
            return res
            
    def run(self,target,references,tolerance=1e-4):
        target_df = self.match_object(target,self.datatable).to_pandas()
        references_df = pd.DataFrame()
        for ref in references:
            ref_df = self.match_object(ref,self.datatable,tolerance=tolerance)
            ref_df['obsmag'] = ref['obsmag']
            references_df = pd.concat([references_df,ref_df.to_pandas()],ignore_index=True,sort=False)

        maglist = []
        for i, ref in references_df.iterrows():
            magi = target_df.MAG_BEST - ref.MAG_BEST + ref.obsmag
            maglist.append(magi)
        mag = np.average(maglist)
        magerr_zp = np.std(maglist) if i>0 else None
        if magerr_zp is None:
            magerr = target_df.MAGERR_BEST.values[0]
        else:
            magerr = np.sqrt(magerr_zp**2 + target_df.MAGERR_BEST.values[0]**2)
        
        return {'mag':mag,'magerr':magerr}
        