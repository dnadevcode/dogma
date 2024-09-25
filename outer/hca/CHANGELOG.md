# Future (4.3):
- Reduce dependence of lldev code (files in bin folder)
implement small things from Wishlist for HCA 4.0 into HCA 4.1, in particular:
- make a script for p-value generation visible, i.e. optional step after computing CC values, assuming that 
p-value database matrix is pre-calculated (otherwise it would take too much time and wouldn't make sense)
- extra functions for additional HCA plots (or ressurect old scripts written for this)
Also todo:
- include UCR_DTW method as optional similarity search method, as well as MP/MPdist (Matrix profile). Note though that 
these two methods are slower (matrix profile runs in O(n^2) as opposed to current method which is O(nlogn), but provides much 
more informative output.
- include possibility to run some sort of SV detection of experiment on the theory (most likely using MP)
- both HCA_GUI and HCA_theory in principle can be run without any gui's now, and all the necessary settings and files can be provided via
'hca_settings.txt', and 'theory_settings.txt'. Need to add clear instructions in the README on how to do this.
- further investigate how parallelization can be employed in the better way. I.e. MASS_PCC is ridiculously easy to parallelize, so one 
option for linear theories would be to concatenate all the theories into one, and then run a single MASS_PCC routine (the indices in each 
subrange of pixels can be identified easily since we know indices of each theory).
- remove some files that are no longer used, i.e. old way to compute theory, unnecessary GUI's, 
- reduce the dependencies, i.e. bin file, because now i.e. the whole +CBT/+Consensus folder is copied over.
- locate some missing files, i.e. gen_dtw_mean for an additional method of kymo alignment, which is currently in a different project, 
/ clustering/src/gen_dtw_mean.m. (why?)

also (copied over from version file)
- will not be merged into gitlab.
- dendogram for consensus
- add additional tab to change/view settings
- change the way stretching is done, so that the end points are not fixed
- add dtw similarity tool
- change that the fft would be computed for smaller series and then stiched together
to avoid memory issues.
- check that consensus is calculated correctly (i.e. zscored!) maybe use matlab's linkage tools
- also add consensus for different kind of maps?


# 5.0.0
- Updated GUI.
#4.6.0 
- Add possibility of importing multi-tif kymos (where frame 2 - kymo, frame 3 -bitmask
- Added (still being updated) spalign algorithm (align method 3), should work better for kymos with parts of edge missing.

#4.5.0
- Added option for minimum length of barcode
- Fix an issue of wrong PCCs when theory barcode is completely flat


# 4.4.1 
small fixes from generating circular plot and filtering

# 4.4
Fixed some issues with MP, added unit tests to see if certain functions are working

# 4.2
- added possibility to generate theory directly via HCA_Gui (if fasta is selected, the program interprets that as needing to generate new theory)
- added compatibility for mp (also with seemingly less bugs when plotting, which appeared before from the way position from the mp is reported)
- added simple p-value generation (random fragments vs non-random theory) in HCA_Run (though it's not pregenerated, so slow for long theories
# 4.1
- Add option to split all barcodes into subbarcodes
- Added extra comparison methods (MP, Spearman), that can be run

# 4.0.2
- Fix theory loading, so that unnecessary folders would not be added to path

# 4.0.0
* Theory generation update
- Change the theory generation function so that it would be more memory efficient. 
- Save theories also as txt files - in the future can keep these and mat file with settings.
- Move the random number generator outside the subfragments generation, so that we wouldn't have
repetitive structures. So theories with NNNNN's will look slightly different
- Added extra parameters for theory generation (i.e. overlap length, how many fragments to cut the theory into)
which were previously hardcoded but now can be changed through the settings file
- Removed output GUI, so all the necessary outputs are saved to user provided location on disk
* Comparison of theory to experiments
- Included a more memory-efficient code, MASS_PCC, as additional possibility to speed up calculations (especially against long theories)
This can be enabled using 
comparisonMethod = 'mass_pcc' in the files/hca_settings.txt; % method for comparison of exp vs theory

- re-enabled parpool over the theories (so should ideally give speed up equal to number of processors (of course in reality less 
because of overhead etc)
- 

# 3.2.8.4

- add a check whether there is any points in the background, if none, we just print nan's as the background values

# 3.2.8.3

- add an option to skip edge detection (for runs where barcode is saved as kymograph/txt file)

# Current version

- Add output table for information scores, and include two options for possible scores
- Change the compare theory to exp function to only allow theory > exp and when exp longer than 20px
- Add subscripts to names of cutout barcodes that have already been defined.
