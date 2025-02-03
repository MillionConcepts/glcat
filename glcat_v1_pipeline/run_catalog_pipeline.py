from gPhoton.pipeline import execute_pipeline
from pyarrow import parquet

metadata = parquet.read_table('/Users/cm/github/glcat/glcat_v1_pipeline/visit_metadata.parquet').to_pandas()

output_dir = "/Users/cm/github/gPhoton2/test_data/"
depth = 120
eclipse = 23456
band = 'NUV'
execute_pipeline(eclipse,band,local_root=output_dir)#,depth=120)
                 
# for eclipse in [36243, 15662, 19197,  4230, 16366, 34213, 23586, 13321,  8614,
#        44175, 15797, 32851, 25981,  9225, 28481, 37312, 18854,  7327,
#         3020, 10740, 40003, 18568, 32289, 40360,  9836, 32941, 35034,
#        20449, 35445, 40057, 41632, 22842, 22962, 32916, 18623, 40408,
#        16409, 18002, 34563, 15231, 12090, 10709, 32259, 16335, 35273,
#        42830, 13757,  7756, 34915,  9181, 13750, 33775, 13856,  4313,
#        14670, 26895, 38910, 13106, 25683,  9261, 28727, 34825, 35407,
#        10509,  2591,  5420, 34522, 13990, 42844, 40402, 17071, 25151,
#        36215, 34126, 22211,  4762, 13424, 13279, 15865, 14679, 25739,
#        13223, 18475, 16355,   931,   842, 40611, 21010, 16258, 34509,
#         5847, 27019,  4246, 13295, 23004,  6792, 32241, 42871,  7853,
#          828]:
for eclipse in [23456,]:
    for mode in [["NUV","mon",None],
                 ["FUV","mof",None],
                 ["FUV","mon",
                  f'{output_dir}e{eclipse}/e{eclipse}-nd-f{str(depth).zfill(4)}-b00-movie-photom-12_8-mon.csv'],
                 ["NUV","mof",
                  f"{output_dir}e{eclipse}/e{eclipse}-fd-f{str(depth).zfill(4)}-b00-movie-photom-12_8-mof.csv"]]:
        #Run NUV photometry with extraction on NUV
        execute_pipeline(
            eclipse,
            mode[0],
            depth=depth,
            threads=4,
            local_root=output_dir,
            recreate=False,
            aperture_sizes=[1.5, 2.3, 3.8, 6.0, 9.0, 12.8, 17.3, 30., 60., 90.],
            write={"movie": True, "image": True},
            coregister_lightcurves=True,
            photometry_only=False,
            compression="rice",
            suffix=mode[1],
            source_catalog_file=mode[2],
        )

#from glcat_merge_utils import *
#from lightcurve_utils import *
from glcat_catalog_utils import accumulate_run_data, make_catalog
from glcat_merge_utils import crossmatch_catalogs
from pyarrow import parquet
# note that this is a different metadata file from the one above!
metadata = parquet.read_table('/Users/cm/github/gPhoton2/gPhoton/aspect/metadata.parquet').to_pandas()

data = accumulate_run_data(23456,metadata=metadata)
ncat = make_catalog(data,'NUV')
fcat = make_catalog(data,'FUV')
cats = {'NUV':ncat,'FUV':fcat}
