# EEG-power-spectrum
Compute the power spectrum of an EEG file.

The power spectrum of an EEG file can be computed using the `tf_ofEEG` function.

This function returns the frequency power of a provided EEG file, 
computing using a tapered sliding time window. The power spectrum 
of each time window is then averaged across time to identify the 
overall power spectrum for the entire time course. 

This function also includes basic outlier detection in order to 
identify outlier timepoints to be removed prior to averaging 
across the timecourse. This is done by transforming each frequency 
band into z-scores, computed using either standard deviation ('sd')
or median absolute deviation ('mad') (see 'ztrans_vector.m'), and 
flagging any timepoints within each frequency band which exceed a 
user-specified outlier threshold (e.g., zthresh = 3). Any timepoints 
whose absolute z-score exceeds this threshold is then flagged as an 
outlier within its frequency band. The percentage of frequency bands 
in which timepoints are flagged as outliers is then computed. If the 
proportion of frequency bands on which a timepoint is flagged exceeds 
the user-specified threshold (e.g., pthresh = .95), then the timepoint 
is marked for removal. 

Please note that such outlier detection methods are very simple, and 
should not be used as the sole means out outlier detection. 

  `eeg_file`      = 'Filepath for the EEG file to be analyzed. Must be .mat format.'
  
  `cfg`           = configuration parameters to be used in the analysis
  
  `cfg.tapper`    = The taper to be used for the sliding time window ('hann', 'hamming', 'gauss')
  
  `cfg.timewin`   = The width of the sliding time window, in ms.
  
  `cfg.winstep`   = How much the sliding time window should move across the EEG time course with each iteration, in ms.
  
  `cfg.zthresh`   = The z-score threshold to be used in outlier detection. Must be a positive number.
  
  `cfg.pthresh`   = The proportion of frequency bands on which a timepoint must exceed the z-threshold to be marked as an outlier for removal. Must be a number between 0 and 1.
  
  `cfg.outmethod` = The deviation statistic to be used when computing z-scores for outlier detection ('sd' or 'mad').
  
  `cfg.freq_chans` = An array containing the channels to plot. 
