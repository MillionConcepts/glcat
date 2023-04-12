import subprocess
import pandas as pd
import io
from util import make_refined_aspect_table, make_file_names
from typing import Optional


def execute_combine(
                    eclipse,
                    expt,
                    num_frames,
                    run_type,
                    dose: Optional[bool] = False,
                    threshold: Optional[float] = None,
                    star_size: Optional[int] = None,
                    output_directory: Optional[str] = "/home/bekah/glcat_tests/"
                    ):
    opt = {
        'eclipse': eclipse,
        'expt': expt,
        'num_frames': num_frames,
        'threshold': threshold,
        'star_size': star_size,
        'dose': dose,
        'run_type': run_type,
        'output_dir': output_directory
    }
    file_names = make_file_names(opt)
    asp_table = make_refined_aspect_table(file_names)
    asp_table.to_csv(file_names["asp_df"])
    return
