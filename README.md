# Streak-data-plotting
Plotting and corrections for streak data from the PK Lab at McGill.<br/>

For the streak data GUI see: https://github.com/dstrande/sci-streak

Only tested on Windows with Python 3.10.

Requires:  numpy, pandas, scipy, matplotlib, h5py
Recommended: lmfit, numdifftools, jupyterlab

Using "conda env create -f streak.yml" with the included yaml file will install most of the depedencies.

scikit-image has a useful function called block_reduce. This can also be done manually with a custom function.

<br/>

Scripts:

  - The experiment.py file should be used first to fill the experiment.json file. This contains all the experimental parameters for each streak file.
  - bksub_h5_ps.py is used to create the hdf5 file from the raw txt data and background files (for picosecond data).
  - bksub_h5_ns.py is used to create the hdf5 file from the raw txt data and background files (for nanosecond data).
  - streak.py creates 3 simple plots: full energy and time traces, full spectra and full decays. 
  - spectra(t).py creates spectra at varying times of the trace (not currently included, needs to be updated a lot).
  - corr_spec and exp_irf_fit.py are modules used to correction spectra and to fit multiexponential equations respectively. 

<br/>

Folders:

  - data: where all data such as extracted fwhm txt files and the .hdf5 file go.
  - decays: where decay plots go.
  - spectra: where spectra plots go.
  - scripts: where any scripts to run go.
  - spectra(t): where any spectra plot along time go.
  - spectra_fits: where any fits of spectra along time go.
  - decay_fits: where any fits of the decays go.

<br/>
