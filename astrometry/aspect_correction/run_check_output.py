from combine_wcs import execute_combine

if __name__ == "__main__":
    execute_combine(9869, 1, 1639, 'xylist', dose=True, threshold=0.8, star_size=3)