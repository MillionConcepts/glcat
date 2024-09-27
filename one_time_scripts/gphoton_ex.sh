#!/bin/bash

conda activate gphoton2

python pipeline_cli.py {eclipse} {band}" \
               f" --threads=4 --local_root={local_root} --stop_after='photonpipe'" \
               f" --compression=rice --write={{'movie':False,'image':False}} --extended_photonlist=True --aspect='aspect' "
