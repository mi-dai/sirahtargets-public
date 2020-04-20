import pandas as pd
def get_photometry_from_web(data,datamjd):
    filters = [x.split('"',maxsplit=1)[1].split('"')[0] for x in data.split('label')[1:]]
    phot = [x.split('[[',maxsplit=1)[1].split(']]')[0].replace(']','').split('[') for x in data.split('data')[1:]]
    df_list = []
    for f,p in zip(filters,phot):
        for arr in p:
            df = {}
            time,mag,magerr = arr.strip().split(',')[:3]
            mjd = float(time)+datamjd
            df['mjd'] = mjd
            df['mag'] = float(mag)
            df['magerr'] = float(magerr)
            df['filter'] = f.split()[-1][0]
            df_list.append(df)
    return pd.DataFrame(df_list)