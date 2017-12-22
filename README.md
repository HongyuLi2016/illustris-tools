# illustris-tools
Utilities for downloading and analyzing the illustris cosmological simulation data.  
See [Illustris website](http://illustris-project.org/) for more infomation about the simulation and the data format.
Please Note that all the code only works on __python2__.
### Contents
* `data`: download cutout file (in hdf5 file format) for subhalo, halo and merger tree.
  * `illustris-get_cutout.py`: __python illustris-get_cutout.py -s 135 -i 3__ will download the particle cutout file of subhalo3 at sanpshot 135 in a folder named snap135-subhalo3 (only some necessary data field of hdf5 file will be donwloaded. Please use illustris-get_cutout_all.py if you need all the data field.)
  * `illustris-get_tree.py`: __python illustris-get_tree.py subhalo3__ will download the merger tree of subhalo3 in a folder named subhalo3.
* `mock`: make mock images, MGEs and Voronoi binned IFU data for a subhalo. __MGE_FIT_SECTORS__ and __Voronoi binning__ packages are required (available here http://www-astro.physics.ox.ac.uk/~mxc/software/). 
  * `illustris-make_img.py`
  * `illustris-make_mge.py`
  * `illustris-make_ifu_bins.py`
  * `illustris-make_ifu_data.py`
  * `illustris-dump_ifu_data.py`
  * `illustris-v_plot.py`
* `shape`: measure shape at different radii for a subhalo.
* `tree`: analyze merger tree of a subhalo and plot some basic figures.
* `utils`: surpport functions for the modelus listed above.
