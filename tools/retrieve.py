import subprocess


def get_eclipse_data(eclipse):
    """ for running on ec2 w s3 bucket backplanes mounted """
    eclipse_pad = "e"+str(eclipse).zfill(5)
    #TODO: set up a folder for processing that's not test_data
    cmd = f"aws s3 cp s3://backplane_test/{eclipse_pad} gPhoton2/test_data --recursive"
    subprocess.call(cmd, shell=True)
    return