# Streak-data-plotting
Plotting and corrections for streak data from the PK Lab at McGill.<br/> Does not include corrections to time and wavelength necesasry for raw files.  

Only tested on Windows.

Requires:  numpy, pandas, scipy, matplotlib, scikit-image, numdifftools, h5py

<br/>

Scripts:

  - make_h5.py is used first to create the hdf5 file from txt data. 

  - streak.py creates 3 simple plots: full energy and time traces, full spectra and full decays. 

  - spectra(t).py creates spectra at varying times of the trace.

  - spectra(t)\_fit.py fits spectra along time with a single gaussian.

  - spectra(t)\_fit_multiple.py fits multiple gaussians to spectra along time. 

<br/>

Files you have to add:

  - experiment.txt contains information about the experiment.

  - names.txt is a file containing information about each streak dataset (i.e. excitation wavelength and fluence).

