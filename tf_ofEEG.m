function [freq_pwr] = tf_ofEEG(eeg_file, cfg)
%myFun - Compute the frequency power of a provided EEG .mat file.
%
% Syntax: [TFRhann] = tf_ofEEG((eeg_file, cfg_specifications)
%
% This function returns the frequency power of a provided EEG file, 
% computing using a tapered sliding time window. The power spectrum 
% of each time window is then averaged across time to identify the 
% overall power spectrum for the entire time course. 
% 
% This function also includes basic outlier detection in order to 
% identify outlier timepoints to be removed prior to averaging 
% across the timecourse. This is done by transforming each frequency 
% band into z-scores, computed using either standard deviation ('sd')
% or median absolute deviation ('mad') (see 'ztrans_vector.m'), and 
% flagging any timepoints within each frequency band which exceed a 
% user-specified outlier threshold (e.g., zthresh = 3). Any timepoints 
% whose absolute z-score exceeds this threshold is then flagged as an 
% outlier within its frequency band. The percentage of frequency bands 
% in which timepoints are flagged as outliers is then computed. If the 
% proportion of frequency bands on which a timepoint is flagged exceeds 
% the user-specified threshold (e.g., pthresh = .95), then the timepoint 
% is marked for removal. 
% 
% Please note that such outlier detection methods are very simple, and 
% should not be used as the sole means out outlier detection. 
%
%   eeg_file = 'Filepath for the EEG file to be analyzed. Must be .mat format.'
%   cfg      = configuration parameters to be used in the analysis
%   cfg.tapper = The taper to be used for the sliding time window ('hann', 'hamming', 'gauss')
%   cfg.timewin = The width of the sliding time window, in ms.
%   cfg.winstep = How much the sliding time window should move across the EEG time course with each iteration, in ms.
%   cfg.zthresh = The z-score threshold to be used in outlier detection. Must be a positive number.
%   cfg.pthresh = The proportion of frequency bands on which a timepoint must exceed the z-threshold to be marked as an outlier for removal. Must be a number between 0 and 1.
%   cfg.outmethod = The deviation statistic to be used when computing z-scores for outlier detection ('sd' or 'mad').
%   cfg.freq_chans = An array containing the channels to plot. 
% 
% See also: ztrans_vector.m and outlier_tpnts

clc

%% test variables
cd('C:\Users\enter\OneDrive - University of Toronto\Labs\Inzlicht Lab\outside_project\data\lvl1\RawData\1030\eegToMatPreprocessed_ICA_cleaned');
eeg_file = 'john_audio_1030_sr_continuous_icacleaned.mat';

cfg.tapper  = 'hann';
cfg.timewin = 500; % in ms
cfg.winstep = 100; % in ms

cfg.zthresh   = 3;   % for outlier timepoint detection
cfg.pthresh   = .95; % for outlier timepoint detection
cfg.outmethod = 'mad';

cfg.freq_chans = {'Pz'}; % channels to plot

%% load in data for analysis
load(eeg_file);
[filepath, filename, fileext] = fileparts(eeg_file);

% channel index
cfg.chanlocs = EEG.chanlocs;

%% transform time to idx
timewinidx = round(cfg.timewin/(1000/EEG.srate));
winstepidx = round(cfg.winstep/(1000/EEG.srate));
pretrialt  = -50; % in ms

%% default variables
tapper.hann    = .5*(1-cos(2*pi*(0:timewinidx-1)/(timewinidx-1)));
tapper.hamming = .54 - .46*cos(2*pi*(0:timewinidx-1)/(timewinidx-1));
tapper.gauss   = exp(-.5*(2.5*(-timewinidx/2:timewinidx/2-1)/(timewinidx/2)).^2);

%% linear detrend data
d = detrend(EEG.data'); % TODO: ask Hause is should be applied to all channels

%% begin fft of time windows
[~,stime] = min(abs(EEG.times - pretrialt)); % start time
stime = stime + timewinidx/2;
ntimewin  = floor((EEG.pnts - timewinidx - stime/2) / winstepidx); % number of time windows

% frequency list of spectral power plot
f = linspace(0,EEG.srate/2,floor(length(tapper.hann)/2)+1); % frequencies of FFT

% loop through each channel
for chan = 1%:EEG.nbchan

    % create TF matrix 
    tf = zeros(floor(length(tapper.hann)/2),EEG.pnts, 1);

    % loop through time points
    for twin = 1:ntimewin

        if twin ~= 1
            stime = stime + winstepidx;
        end

        dfft = fft(d(stime:stime+timewinidx-1,chan)' .* tapper.hann); % fast fourier transformation of the detrended time series within a timewin

        % TF matrix input column of data at selected time point
        selected_time = stime-timewinidx/2:stime+timewinidx/2;
        tf(:,selected_time,chan) = repmat(abs(dfft(2:floor(length(tapper.hann)/2)+1)')*2,1,length(selected_time));
    end

    % columns with all zeros to remove
    cols_with_all_zeros       = find(all(tf==0));
    tf(:,cols_with_all_zeros) = [];
    
    % detect outliers
    outliers       = outlier_tpnts(tf, cfg.zthresh, cfg.pthresh, cfg.outmethod);
    tf(:,outliers) = [];

    % remove deleted times
    freq_times                                  = EEG.times;
    freq_times([cols_with_all_zeros, outliers]) = [];

    figure
    imagesc(freq_times,f,log10(tf+1))
    set(gca,'clim',[-1 1]*4)

    freq_pwr = mean(tf,2); % spectral power
    figure
    plot(freq_pwr)
end

end