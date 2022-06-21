# Streak-data-plotting
Plotting and corrections for streak data from the PK Lab at McGill.<br/>

Only tested on Windows.

Requires:  numpy, pandas, scipy, matplotlib, scikit-image, h5py
Recommended: lmfit, numdifftools

<br/>

Scripts:

  - The experiment.py file should be used first to fill the experiment.json file. This contains all the experimental parameters for each streak file.
  - bksub_h5.py is used to create the hdf5 file from the raw txt data and background files. 
  - streak.py creates 3 simple plots: full energy and time traces, full spectra and full decays. 
  - spectra(t).py creates spectra at varying times of the trace.
  - corr_spec and exp_irf_fit.py are modules used to correction spectra and to fit multiexponential equations respectively. 

<br/>

Folders:

  - data: where all data such as extracted fwhm files and the .h5 go.
  - decays: where decay plots go.
  - spectra: where spectra plots go.
  - scripts: where any scripts to run go.
  - spectra(t): where any spectra plot along time go.
  - spectra_fits: where any fits of spectra along time go.
  - decay_fits: where any fits of the decays go.

<br/>
