""" this was in a branch of gphoton, no longer necessary but may be useful someday """

def add_slewframes(eclipse, band, ):
    """adds slewframes to aspect table for eclipse so that the slew
    frames can be added to photonlists"""
    from pyarrow import parquet
    from gPhoton.io.mast import get_raw_paths, download_data
    from astropy.io import fits

    # loading aspect table to add time stamp and flags
    parq = parquet.read_table('/home/bekah/gphoton_working/gPhoton/aspect/aspect.parquet')
    aspect = parq.to_pandas()
    aspect = aspect[aspect["eclipse"] == eclipse]
    aspect = aspect.reset_index()

    # download scst
    scst = fits.open(download_data(eclipse, "scst", band))
    scst = fits.open(download_data(eclipse, "scst", band))
    scst_pd = pd.DataFrame(scst[1].data)
    scst_pd = scst_pd.rename(columns={"pktime": "time"})

    # merge to get slew frames
    slew_frames = scst_pd.merge(aspect, how='left', on=['time'])
    slew_frames = slew_frames[slew_frames['hvnom_nuv'] == 1]

    # use acs solns for ra / dec / roll for slew frames in asp soln
    slew_frames['ra'] = slew_frames['ra'].fillna(slew_frames['ra_acs'])
    slew_frames['dec'] = slew_frames['dec'].fillna(slew_frames['dec_acs'])
    slew_frames['roll'] = slew_frames['roll'].fillna(slew_frames['roll_acs'])

    slew_frames = slew_frames.rename(columns={"eclipse_x": "eclipse", "flags_y": "flags"})
    slew_frames = slew_frames.filter(['eclipse', 'time', 'flags', 'ra', 'dec', 'roll'],
                               axis=1)
    slew_frames = slew_frames.astype({'ra': 'float64',
                                      'dec': 'float64',
                                      'roll': 'float64'})
    slew_frames = slew_frames.reindex(range(len(slew_frames)))

    # save to parquet
    #TODO: change dir name here
    slew_frames.to_parquet('/home/bekah/gphoton_working/gPhoton/aspect/aspect2.parquet', compression=None)

    return
