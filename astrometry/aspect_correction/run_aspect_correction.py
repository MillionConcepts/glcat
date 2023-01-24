from main import execute_refiner

if __name__ == "__main__":
    execute_refiner(
        9869, 1, 1639, 'xylist', dose=True, threshold=0.8, star_size=3
    )
