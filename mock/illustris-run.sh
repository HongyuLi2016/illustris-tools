#!/bin/bash
# File              : illustris-run.sh
# Author            : Hongyu Li <lhy88562189@gmail.com>
# Date              : 22.12.2017
# Last Modified Date: 22.12.2017
# Last Modified By  : Hongyu Li <lhy88562189@gmail.com>

# make mock image -p phi -i inclination -o oblate rotatorion
illustris-make_img.py $1 -p $2 -i $3 -o  
# make stellar mge using the mock image created by illustris-make_img.py
illustris-make_mge.py $1
# Voronoi bin
illustris-make_ifu_bins.py $1
# Calculate Vel, Vdisp and their error in each Voronoi bin
illustris-make_ifu_data.py $1
# dump IFU data to a fits file and plot velocity map
illustris-dump_ifu_data.py $1
illustris-v_plot.py $1
