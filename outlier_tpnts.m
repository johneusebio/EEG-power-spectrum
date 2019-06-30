function outliers = outlier_tpnts(tf, zthresh, pthresh, method)
%outlier_tpnts - detect outlier timepoints in a EEG time-frequency matrix
%
% Syntax: outliers = outlier_tpnts(tf)
%
% This function detects outlier timepoints in a EEG time-frequency matrix
% and returns the column indices.
% 
% This is accomplished by computing teh z score of each time point within 
% all frequency bands (i.e., rows) and finding time-points which exceed 
% the user-specified z-threshold. 
% 
% The time points which are found to be outliers more frequently than the 
% user-specified probability threshold will be marked as outliers.

zmat   = arrayfun(@(ROWIDX) ztrans_vector(tf(ROWIDX,:), method), (1:size(tf,1)).', 'UniformOutput', false);
% disp(zmat);
zmat   = cell2mat(zmat);
outmat = abs(zmat) > zthresh;
pmat   = mean(outmat, 1);
    
outliers = pmat > pthresh;
outliers = find(outliers);

end