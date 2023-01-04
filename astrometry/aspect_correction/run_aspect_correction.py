from main import execute_refiner

if __name__ == "__main__":
    execute_refiner(
        21442, 1, 161, 'xylist', threshold=.8, star_size=2, verify=True
    )
