from new_aspect_pipeline import refine_eclipse

refine_eclipse(eclipse=1532,
               xylist=False,
               band="NUV",
               backplane_root="/media/bekah/BekahA/backplanes",
               aspect_root="/media/bekah/BekahA/aspect",
               gphoton_root="/home/bekah/gPhoton2",
               ext="gzip",
               threads=4)


