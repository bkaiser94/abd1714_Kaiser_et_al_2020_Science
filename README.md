# abd1714_Kaiser_et_al_2020_Science
Reduction and Analysis Code for Kaiser et al. 2020, Exo-Planetesimal Lithium Pollution of a White Dwarf

To run this code you have to retrieve the thin-H C/O-core WD cooling model tables from the Montreal website LINK TO BE INSERTED, and then compile them into a single CSV file with the corresponding name in the scripts.

You'll also need to retrieve the relevant table of A(Li) for the F and G stars from Bensby and Lind 2018 and the relevant table from Bensby, Feltzing, and Oey 2014 for the other abundances and ages. These two tables will then need to be merged into a single CSV file with the name from the scripts.

You'll need to get the model atmospheres from Moehler et al. 2014 as well for the standards for flux-calibration.

You'll also need to get the extinction curve data from Stritzinger et al. 2005.


You'll also have to change a lot of paths that are hard-coded in various places to my own computer (e.g. /BenKaiser/Desktop/) obviously won't work on your computer. Basically, you're gonna be sitting there for awhile pulling your hair out if you want to run this locally. Just like I had to do to make it.

All scripts *should* be compatible with Python 3, but I personally usually run the reduction scripts in Python 2.7.15. I make no guarantees on compatibility with various versions of packages as well.

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