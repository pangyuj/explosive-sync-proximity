% This function simulates coupled Stuart-Landau model with time delay,
% coherent external stimulation and Adaptive Feedback, Improved Euler method
%
% Inputs
% MAT: Connectivity matrix
% dist: Distance matrix
% speed: Speed of interaction (delay is defined by dist/speed)
% noise: Gaussian noise for each step
% W: initial frequencies
% T: time vector for simulation
% zeta: parameter for adaptive feedback. zeta=0 means no feedback effect 
% and higher value (zeta=1,2,3,4 ...) means stronger effect.
% initV: initial conditions (initV = 2*pi*rand(1,sizeMAT*nE)-pi)
% u: external perturbation to increase Z along with real axis. same length with T
%   'dist', 'speed', 'W' and 'T' should be in same scale (ex: meter, meter/sec, Hz, sec)
%
% Outputs
% T: time corresponds to XO or XE
% ZT: signal time series
%
% 2016.11.30. modified by Hyoungkyu Kim to improve calculation speed
% Last updated for gpu by Joseph Lee 12/1/2016
% updated by JY Moon 3/11/2017
% Reorganized with ImprovedEuler and added annotations by Pangyu Joo 9/16/2024

function [T,ZT] = IE_stuartlandau_distdelay_stim_af(MAT, dist, speed, noise, W, T, zeta, initV ,u)
alpha = 1; % parameter for the Stuart-Landau oscillator
%% parameters for simulation
if ~exist('MAT','var')
    MAT = 0; % adjacency matrix of the network
end
sizeMAT=length(MAT); % number of nodes in the network
if ~exist('dist','var')
    dist = zeros(sizeMAT); % distance matrix of the network
end
if ~exist('speed','var')
    speed=5; % speed of signal propagation
end
if ~exist('noise','var')
    noise=0; % noise level
end
if ~exist('W','var')
    w0=10; % base frequency
    W = w0+randn(1,sizeMAT); % frequencies of the oscillators
end
if ~exist('T','var')
    T = 0:0.001:10; % time vector
end
if ~exist('zeta','var')
    zeta = 0; % parameter for adaptive feedback
end
if ~exist('initV','var')
    initV = randn(1,sizeMAT)+1i*randn(1,sizeMAT); % initial conditions
end
if ~exist('u','var')
    u= zeros(size(T)); % external stimulation
end

%% initialization
W = W(:).*2*pi; % convert frequencies to angular frequencies
N = length(T); % number of time steps

[m_x, m_y, Coupling_str] = find(MAT); % find non-zero components of the adjacency matrix
idx=find(MAT); % linear indices for non-zero components

dt=(T(end)-T(1))/(length(T)-1);
delay = dist./speed; % delay times due to signal propagation
delay = round(delay/dt); % round delay times to nearest integer
delaymax=max(max(delay)); % maximum delay time
ZT = zeros(N,sizeMAT); % preallocate memory for the solution
ZT(1,:) = initV; % set initial conditions
Z_temp = ZT(1,:).';
delay_time=min(delaymax,N-1); % maximum delay time
noise_series=noise.*randn(sizeMAT,2); % assume noise is a function of time (same noise value for same time)
%% simulation before maximum delay time (no delay)
for ti=1:delay_time
    % constant parameters
    h=T(ti+1)-T(ti);
    SumConn=0; % interaction term is set to 0 before maximum delay time
    u_in=0; % perturbation term is set to 0 before maximum delay time

    % Improved Euler method
    k1=h*EOM_SL(Z_temp, alpha, W, SumConn, noise_series(:,1), u_in, zeta);
    k2=h*EOM_SL(Z_temp+k1, alpha, W, SumConn, noise_series(:,2), u_in, zeta);

    Z_temp=Z_temp+(k1+k2)/2;
    ZT(ti+1,:)=Z_temp.'; % update ZT

    noise_series(:,1)=noise_series(:,2);
    noise_series(:,2)=noise.*randn(sizeMAT,1);
end
%% simulation after maximum delay time (with delay)
ZTj_term2=N*repmat([0:1:sizeMAT-1]',1,sizeMAT);
ZTj_term2=ZTj_term2(idx);  %linear indices for connected couplings

delay2=delay(idx); % delays for connected couplings

for ti=delaymax+1:N-1
    h=T(ti+1)-T(ti);
    % Calculate delayed terms for the differential equations (for k1)
    ZT_list = Coupling_str.*ZT(ti-delay2+ZTj_term2);
    ZT_mat = sparse(m_x,m_y,ZT_list,sizeMAT,sizeMAT);
    SumConn=sum(ZT_mat).'; %delayed interaction term (summed)
    k1=h*EOM_SL(Z_temp, alpha, W, SumConn, noise_series(:,1), u(ti), zeta);
    
    %update ZT with k1
    ZT(ti+1,:)=(Z_temp+k1).'; 

    % Calculate delayed terms for the differential equations (for k2) 
    ZT_list = Coupling_str.*ZT(ti+1-delay2+ZTj_term2);
    ZT_mat = sparse(m_x,m_y,ZT_list,sizeMAT,sizeMAT);
    SumConn=sum(ZT_mat).'; %delayed interaction term (summed)
    k2=h*EOM_SL(ZT(ti+1,:).', alpha, W, SumConn, noise_series(:,2), u(ti+1), zeta);

    %improved Euler
    Z_temp=Z_temp+(k1+k2)/2;
    ZT(ti+1,:)=Z_temp.'; % update ZT with final result

    noise_series(:,1)=noise_series(:,2);
    noise_series(:,2)=noise.*randn(sizeMAT,1);
end

end
%% functions
% This function defines the differential equations for the Stuart-Landau oscillators.
function dx = EOM_SL(Z, alpha, W, SumConn, noise_in, u, zeta)
% adaptive feedback
[Ra] = LocalAllOrderParameter(Z);
RaZ = Ra.^(zeta);
% equation of motion
dx = (alpha+1i*W-Z.*conj(Z)).*Z + RaZ.*SumConn + u + noise_in;
end

% This function calculates Local synchrony of each nodes from complex time series.
% 10.31.2016 Joon-Young Moon
% modified to receive complex input by Pangyu Joo (9.16.2024)
function [Ra] = LocalAllOrderParameter(Z)
Z=Z./abs(Z);
mean_i = mean(Z,'omitnan');
mean_i = (mean_i + Z )/2;
Ra= abs(mean_i);
end