#!/bin/bash
illustris-make_img.py $1 -i 90
illustris-make_mge.py $1
illustris-make_ifu_bins.py $1
illustris-make_ifu_data.py $1
illustris-dump_ifu_data.py $1
illustris-v_plot.py $1
