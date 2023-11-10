from test_aspect_soln import run_compare

eclipses = [580, 1696, 1687, 1710, 1752, 1817, 1777, 6411,
     5267, 7535, 1819, 6268, 6411, 6413, 10507, 11604, 15466, 23725]

for e in eclipses:
    try:
        print(f"Testing eclipse {e}")
        run_compare(e, 0, "NUV", runnote="", local_root="/home/bekah/gPhoton2/")
    finally:
        print("done with test.")