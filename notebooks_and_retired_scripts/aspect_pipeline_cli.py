"""
simple command-line interface to execute_pipeline()
"""
import fire

from notebooks_and_retired_scripts.main_compare_asp import run_compare

# tell fire to handle command line call
if __name__ == "__main__":
    fire.Fire(run_compare)

