from test_aspect_soln import run_compare

if __name__ == "__main__":
    for i in [1687,1777,5267,23725,23686,15466]:
        try:
            run_compare(
            i, 1, "NUV", runnote="", local_root="/home/bekah/gPhoton2/")
        except:
            print("something no work :(")