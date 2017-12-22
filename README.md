# illustris-tools
Utilities for downloading and analyzing the illustris cosmological simulation data.  
See [Illustris website](http://illustris-project.org/) for more infomation about the simulation and the data format.
Please Note that all the code here only works under __python2__.
### Contents
* `data`: download cutout file (in hdf5 file format) for subhalo, halo and merger tree.
  * `illustris-get_cutout.py`: __python illustris-get_cutout.py -s 135 -i 3__ will download the particle cutout file of subhalo3 at sanpshot 135 in a folder named __snap135-subhalo3__ (only some necessary data field of hdf5 file will be donwloaded. Please use illustris-get_cutout_all.py if you need all the data field.)
  * `illustris-get_tree.py`: __python illustris-get_tree.py subhalo3__ will download the merger tree of subhalo3 in a folder named subhalo3.
* `mock`: make mock images, MGEs and Voronoi binned IFU data for a subhalo. _MGE_FIT_SECTORS_ and _Voronoi binning_ packages are required (available here http://www-astro.physics.ox.ac.uk/~mxc/software/). 
  * `illustris-make_img.py`: __python illustris-make_img.py snap135-subhalo3 -p 30 -i 90 -o__ will create mock images for a subhalo (cutout data file need to be downloaded first). -p -i are rotation angle phi and inclination, -o means using oblate rotation. Please read the code or ask the author for more details.  The output files will be in a folder named imgs in snap135-subhalo3.
  * `illustris-make_mge.py`: __python illustris-make_mge.py snap135-subhalo3__ will create MGE using the mock image (img_M.npy). The output files will be in a folder named mge in snap135-subhalo3.
  * `illustris-make_ifu_bins.py`: __python illustris-make_ifu_bins.py snap135-subhalo3__ will create Voronoi bins for the IFU data. The output files will be in snap135-subhalo3/ifu
  * `illustris-make_ifu_data.py`: __python illustris-make_ifu_data.py snap135-subhalo3__ will calculate mean velocity, velocity dispersion and their errors for the particles in each Voronoi bin. The output files will be in snap135-subhalo3/ifu
  * `illustris-dump_ifu_data.py`: __python illustris-dump_ifu_data.py snap135-subhalo3__ will dump the IFU data into a fits file IFU.fits in snap135-subhalo3/ifu
  * `illustris-v_plot.py`: __python illustris-v_plot.py snap135-subhalo3__ will read the IFU.fits file and plot the velocity maps. The output figure will be in snap135-subhalo3/ifu
* `shape`: measure shape at different radii for a subhalo.
* `tree`: analyze merger tree of a subhalo and plot some basic figures.
* `utils`: surpport functions for the modelus listed above.
