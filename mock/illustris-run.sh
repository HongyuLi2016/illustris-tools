#!/bin/bash
illustris-make_img.py $1 -p $2 -i $3 -o
illustris-make_mge.py $1
illustris-make_ifu_bins.py $1
illustris-make_ifu_data.py $1
illustris-dump_ifu_data.py $1
illustris-v_plot.py $1
