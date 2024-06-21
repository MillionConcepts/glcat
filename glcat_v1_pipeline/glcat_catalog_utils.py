import numpy as np
import pandas as pd
from lightcurve_utils import counts2mag

def accumulate_run_data(eclipse,metadata=None,
                        datadir = "/Users/cm/github/gphoton2_refactor/gPhoton2/test_data/",
                        apers = [1.5, 2.3, 3.8, 6.0, 9.0, 12.8, 17.3],
                       ):
    edir = f"e{str(eclipse).zfill(5)}"
    data = {'eclipse':eclipse,}
    if not metadata is None:
        data['obstype']=metadata[metadata['eclipse']==eclipse]['obstype'].values[0]
    at = {}
    for p in 'nf':
        for m in 'nf':
            at[f'{p}{m}'] = {}
            for a in apers:
                at[f'{p}{m}'][a] = pd.read_csv(
                    f"{datadir}{edir}/{edir}-{p}d-b00-f0120-movie-photom-{str(a).replace('.','_')}-mo{m}.csv",
                    index_col=None)
    data['phot'] = at

    data['expt'] = {'NUV':pd.read_csv(f"{datadir}{edir}/{edir}-nd-b00-f0120-movie-exptime.csv",index_col=None),
                    'FUV':pd.read_csv(f"{datadir}{edir}/{edir}-fd-b00-f0120-movie-exptime.csv",index_col=None)}
    return data

def derived_photometry_table(data,aper,code):
    band = {'n':'NUV','f':'FUV'}[code[0]]
    expt = data['expt'][band]['expt'].sum()
    cps = np.array(data['phot'][code][aper]['aperture_sum'].tolist())/expt
    cps_err = np.sqrt(data['phot'][code][aper]['aperture_sum'].tolist())/expt
    mag = counts2mag(cps,band)
    mag_err_upper = np.abs(counts2mag(cps-cps_err,band) - mag)
    mag_err_lower = np.abs(counts2mag(cps+cps_err,band) - mag)

    return pd.DataFrame(
        {'CPS':cps,'CPS_ERR':cps_err,
         'MAG':mag,'MAG_ERR_UPPER':mag_err_upper,'MAG_ERR_LOWER':mag_err_lower
         #'MASK_FLAG':(data['phot'][code][aper]['aperture_sum_mask'] != 0).astype(int),
         #'EDGE_FLAG':(data['phot'][code][aper]['aperture_sum_edge'] != 0).astype(int),
        }
    )

def accumulate_photometry(data,code,
                          apers=[1.5, 2.3, 3.8, 6.0, 9.0, 12.8, 17.3],
                         ):
    band = {'n':'NUV','f':'FUV'}[code[0]]
    phot = pd.DataFrame()
    for i,aper in enumerate(apers):
        suffix = f'_APER{i}'
        phot[f'{band}_SUM{suffix}']=list(data['phot'][code][aper]['aperture_sum'])
        phot[f'{band}_EDGE{suffix}']=list((data['phot'][code][aper]['aperture_sum_edge'] != 0).astype(int))
        phot[f'{band}_MASK{suffix}']=list((data['phot'][code][aper]['aperture_sum_mask'] != 0).astype(int))
        dp = derived_photometry_table(data,aper,code)
        phot = pd.concat([phot,dp.add_prefix(f'{band}_').add_suffix(suffix)],axis=1)
    phot[f'{band}_EXPT']=np.full(len(data['phot'][code][aper]),
                                 data['expt'][band]['expt'].sum())
    return pd.DataFrame(phot)

def make_catalog_by_code(data,code):
    band = {'n':'NUV','f':'FUV'}[code[0]]
    cat = pd.concat([pd.DataFrame(
        {'eclipse':np.full(len(data['phot'][code][12.8]),
                           data['eclipse']),
         'obstype':np.full(len(data['phot'][code][12.8]),
                           data['obstype']),
        }),
                    ],axis=1)
    if len(set(code))==1: # only add eclipse, ra, dec from primary / extracted catalog
        cat = pd.concat([cat,data['phot'][code][12.8][['ra','dec']],
                        ],axis=1)
    cat = pd.concat([cat,
                     data['phot'][code][12.8][['xcenter','ycenter']].add_prefix(f'{band}_'),
                     accumulate_photometry(data,code),
                    ],axis=1)
    if len(set(code))==1: # only the primary / extracted catalog contains extended source information
        cat = pd.concat([cat,pd.DataFrame({f'{band}_EXTENDED':
                                           list(data['phot'][code][12.8]['extended_source'])})],axis=1)
    return cat

def make_catalog(data,band):
    b = band[0].lower()
    code = f'{b}{b}'
    cat_1 = make_catalog_by_code(data,code)
    code = ({"f":"n","n":"f"}[b])+b # there must be a more elegant way to do this
    cat_2 = make_catalog_by_code(data,code)
    cat = pd.concat([cat_1,cat_2],axis=1).rename(columns=str.upper)
    assert cat[['RA']].shape[1]==1
    return cat

