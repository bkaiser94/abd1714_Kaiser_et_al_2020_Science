# abd1714_Kaiser_et_al_2020_Science
Reduction and Analysis Code for Kaiser et al. 2020, Exo-Planetesimal Lithium Pollution of a White Dwarf

This code contains the scripts used to reduce the raw Goodman data to produce the 1-d spectra provided as supplementary data. It also contains the code to produce the figures in the paper (which incorporates some of the analysis). It includes the necessary supporting scripts for the figure production and to produce the input files for the figures (excluding those files that come from other publicly available sources as outlined below). 

The code and reference files refer to "GaiaJ1644-0449" throughout them. This was the previous name we used for WD J1644–0449, and rather than changing it everywhere in the code I just created a function that replaces it in outputs for display purposes.

***************
THINGS YOU WILL HAVE TO DO TO MAKE THIS SCRIPT RUN (THIS IS NOT EXHAUSTIVE IT'S JUST DEFINITE FILES TO RETRIEVE):

To run this code you have to retrieve the thin-H C/O-core WD cooling model tables from the Montreal website http://www.astro.umontreal.ca/~bergeron/CoolingModels/, and then compile them into a single CSV file with the corresponding name in the scripts. This must be done manually; my code does not do this (and did not do it for me).

You'll also need to retrieve the relevant table of A(Li) for the F and G stars from Bensby and Lind 2018 and the relevant table from Bensby, Feltzing, and Oey 2014 for the other abundances and ages. These two tables will then need to be merged into a single CSV file with the name from the scripts (bensby_plotting.py is the script that calls the file). This is another manual task.

You'll need to get the model atmospheres from Moehler et al. 2014 as well for the standards for flux-calibration. This is also a manual task.

You'll also need to get the extinction curve data from Stritzinger et al. 2005. Yet another manual task.

Ah, you'll also need to retrieve the SDSS-BOSS spectral template library from Kesseli et al. 2017 for making the comparison spectra if you want to recreate "telluric corrections" figure from the paper; it's the library behind one of the PyHammer releases. I suppose you could just retrieve the K7 dwarf template as that's the only one used in the Jupyter notebook. You guessed it; must be done manually.

You'll also have to change a lot of paths that are hard-coded in various places to my own computer (e.g. /BenKaiser/Desktop/ obviously won't work on your computer). Basically, you're gonna be sitting there for awhile pulling your hair out if you want to run this locally. Just like I had to do to make it.

All scripts *should* be compatible with Python 3, but I personally usually run the reduction scripts in Python 2.7.15. I make no guarantees on compatibility with various versions of packages as well. Additionally, these scripts were executed at a time and using data obtained during a time that Astropy still had been updating the URL required to get the necessary barycentric corrections; the existing link in the astropy code is now (as of 2020-11-04) no longer able to as accurately deal with barycentric correction stuff for newly obtained data. Again, when the code was executed over a year ago, Astropy did perform corrections the right way. I include this caveat in case anyone tries to use the reduction scripts on their new Goodman data in Python 2.7.


***********

I did not include the files I listed above because they are not mine to distribute but they are very easily accessible to anyone that cares to retrieve them. I hope this also encourages individuals to cite the papers from which the materials originate.


************

The Jupyter notebooks were run in Python 3 with the following package versions (these are not necessarily the required ones; they're just the ones I *know* work because I used them):
-numpy 1.14.3
-matplotlib 2.2.2
-astropy 3.2.3
-scipy 1.3.1

As far as the packages I was running in Python 2.7.15 for the reduction scripts those are as follows:

-numpy 1.11.3
-matplotlib 1.5.3
-astropy 2.0.16
-scipy 1.2.2