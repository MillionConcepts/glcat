from main import execute_refiner

if __name__ == "__main__":
    execute_refiner(
        9869, 32, 50, 'xylist', dose=True, threshold=1.0, star_size=3, verify=False
    )
