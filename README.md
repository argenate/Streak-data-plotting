# Streak-data-plotting
Plotting and corrections for streak data from the PK Lab at McGill.<br/>

Only tested on Windows.

Requires:  numpy, pandas, scipy, matplotlib, scikit-image, h5py
Recommended: lmfit, numdifftools

<br/>

Scripts:

  - The experiment.py file should be used first to fill the experiment.json file. This contains all the experimental parameters for each streak file.
  - backsub_h5.py is used first to create the hdf5 file from the raw txt data and background files. 
  - streak.py creates 3 simple plots: full energy and time traces, full spectra and full decays. 
  - spectra(t).py creates spectra at varying times of the trace.
  - spectra(t)_fit.py fits spectra along time with a single gaussian. For testing purposes and finding fit parameters.
  - spectra(t)_fit_multiple.py fits multiple gaussians to spectra along time. The fit parameters should already be known here.

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
