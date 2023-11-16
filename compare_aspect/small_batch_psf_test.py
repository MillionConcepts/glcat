from test_aspect_soln import run_compare

eclipses = [580]

for e in eclipses:
    try:
        print(f"Testing eclipse {e}")
        run_compare(e, 0, "NUV", runnote="", local_root="/home/bekah/gPhoton2/")
    finally:
        print("done with test.")