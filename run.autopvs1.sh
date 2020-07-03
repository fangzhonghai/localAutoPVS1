#!/bin/bash
in_file=$1
out_file=$2
export PYTHONPATH=/zfssz4/B2C_RD_P2/PMO/fangzhonghai/software/download/:$PYTHONPATH
python3 autopvs1_in_bgianno.py -i $in_file -o $out_file -p 1
