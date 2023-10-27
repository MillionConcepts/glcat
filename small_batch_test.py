from run_correct_on_ec2 import test_eclipse

e = [580, 1436, 1519, 1696, 1687, 1710, 1752, 1817, 1777, 6411,
     5267, 7535, 1819, 6268, 6411, 6413, 10507, 11604, 15466,
     23725, 25191, 25192, 26418, 26420, 1414, 1420, 842, 769,
     589, 7581, 9010]

try:
    for e in list:
        print(f"Testing eclipse {e}")
        test_eclipse(e)
finally:
    print("done with test.")
