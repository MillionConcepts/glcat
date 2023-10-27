from new_aspect_pipeline import refine_eclipse
import subprocess
import os
import shutil
from hostess.profilers import DEFAULT_PROFILER as PRO

eclipse = 815
eclipse_str = 'e' + str(eclipse).zfill(5)

with PRO.context(eclipse_str):
    base = "backplanes"
    base2 = "aspect"
    bucket_path = f"s3://backplanetest/{eclipse_str}/"
    ec2_base = f"/home/ubuntu/{base}"
    ec2_base2 = f"/home/ubuntu/{base2}"
    ec2_full = ec2_base+"/"+eclipse_str

    # get files from s3 bucket
    get_files = f'aws s3 cp "{bucket_path}" "{ec2_full}" --recursive'

    # make dirs on instance
    if not os.path.exists(ec2_base):
        print(f"Creating folder '{base}'...")
        os.makedirs(ec2_base)
    else:
        print(f"Folder '{base}' already exists.")
    if not os.path.exists(ec2_full):
        print(f"Creating subsequent folder '{eclipse_str}'...")
        os.makedirs(ec2_full)
    else:
        print(f"Subsequent folder '{eclipse_str}' already exists.")
    if not os.path.exists(ec2_base2):
        print(f"Creating folder '{base2}'...")
        os.makedirs(ec2_base2)
    else:
        print(f"Folder '{base2}' already exists.")

    try:
        result = subprocess.run(get_files,
                                shell=True,
                                check=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                universal_newlines=True)
        print("Command Output:")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with error: {e}")

    aspect = refine_eclipse(eclipse=eclipse,
                   xylist=False,
                   band="NUV",
                   backplane_root=ec2_base,
                   aspect_root=ec2_base2,
                   gphoton_root="/home/ubuntu/gPhoton2",
                   ext="gzip",
                   threads=4)

    try:
        shutil.rmtree(ec2_base)
        print(f"Deleted folder: {ec2_base}")
    except OSError as e:
        print(f"Error deleting folder: {e}")

print("Pipeline is done running :)")
print(f"It took: {PRO.labels[eclipse_str]['time']} seconds")



