function Z = ztrans_vector(X, method)
%ztrans_vector - z-transform each element of a user-specified vector
%
% Syntax: Z = ztrans_vector(X)
% 
% This function transforms each element of a given numeric vector into
% a z-score. These z-scores are computed using either of two methods:
%   'std' = standard deviation, which is more affected by larger values
%   'mad' = median absolute deviation, which is more robust against such 
%           outliers.

if strcmpi(method, 'std')
    Z = (X - mean(X)) / std(X);
elseif strcmpi(method,'mad')
    Z = (X - mean(X)) / mad(X);
else
    error('The specified outlier detection method must use either ''std'' or ''mad''')
end

end