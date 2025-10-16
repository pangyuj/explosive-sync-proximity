% This function calculates the global order parameter and the order phase from a complex time series.
% Inputs
%   Z: A complex time series (NxT matrix, where N is the number of oscillators and T is the number of time points)
% Outputs
%   or: Global order parameter (average across all time points)
%   or_std: Standard deviation of the order parameter over time
%   or_t: Time-resolved order parameter (order parameter at each time point)
function [or,or_std,or_t] = OrderParameter_Comp(Z)

% Normalize the complex time series to unit vectors (extracting phase information)
phase_vec = Z ./ abs(Z);

% Calculate the mean phase vector across all oscillators at each time point
ph_vector = mean(phase_vec, 2);

% Calculate the time-resolved global order parameter as the magnitude of the mean phase vector
or_vector = abs(ph_vector);

% outputs
or_t = or_vector(1:end); % time-resolved order parameter
or = mean(or_t); % overall global order parameter
or_std = std(or_t); %  the standard deviation of the time-resolved order parameter