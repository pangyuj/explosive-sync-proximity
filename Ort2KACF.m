% This function calculates the kurtosis of the autocorrelation of the
% orderparameter (global synchronization) time series.
% Inputs:
%   or_t: The time-resolved order parameter (vector of synchronization values over time).
%   tau: The lag value for the autocorrelation function (how far back in time to compare values).
%   window_size: The size of the sliding window for calculating autocorrelation (number of time points per window).
%   overlap: overlap of windows (between 0~1)
%
% Output:
%   kacf: Kurtosis of the autocorrelation values over different windows.

function [kacf] = Ort2KACF(or_t, tau, window_size, overlap)

step_size=window_size*(1-overlap);
wn = floor((length(or_t)-window_size) / step_size)+1; % Number of windows
acf_w = zeros(wn, 1); % Initialize autocorrelation array

% Loop through each window and compute autocorrelation at lag 'tau'
for wi = 1:wn
    acf_tmp = autocorr(or_t((wi-1)*step_size + 1 : (wi-1)*step_size + window_size),'NumLags',tau);
    acf_w(wi) = acf_tmp(end); % Store autocorrelation at lag 'tau'
end

% Calculate kurtosis of the autocorrelation values
kacf = kurtosis(acf_w);