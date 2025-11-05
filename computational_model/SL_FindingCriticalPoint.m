% This script simulates a Stuart-Landau network with distance delay and adaptive feedback.
% Calculates the global order parameter, standard deviation, and kurtosis of autocorrelation
% and Pair correlation function (KACF) across different zeta values and coupling strengths,
% to find critical point.
clear; clc; close all;

% set save path and generate folder if not exist
savepath='results_400iter/';
if ~exist(savepath,"dir"); mkdir(savepath); end
shift_peak=1;

% simulation parameters
load('gong78');  % Load connectivity (MAT) and distance matrix (Dmat)
dist = Dmat/1000;  % Convert distance to meters (m)
speed = 7;        % Speed of signal propagation (m/s)
dt = 0.001;       % Time step for the simulation
T = 0:dt:10;     % Time vector (simulation duration 10 seconds to find critical point)
T2 = 0:dt:305;     % Time vector (simulation duration 305 seconds to calculate KACF, KPCF)
Transient_time=5; % Transient time to be removed from analysis
noise = 1;      % Set noise level

% initial conditions
load('W_save');   % Load saved initial frequency values
load('initial_save');  % Load saved initial conditions
ini_ind = 1:400;  % Index for initial conditions used in simulation
W_sets = W_save(ini_ind,:);  % Extract initial frequencies for selected indices
initials = initial_save(ini_ind,:);  % Extract initial conditions for selected indices
tau=0.20; % tau for ACF (in sec)
win_size=10; % window size for ACF and PCF (in sec)
overlap=0.5; % overlap of windows to calculate ACF and PCF [0-1]

% Define zeta values (adaptive feedback strength)
zetas = 0:0.5:3;

% Define coupling strengths
strength = 0:0.05:5;

% Loop over each zeta value (adaptive feedback strength)
fprintf('Finding critical points \n')
for zi = 1:length(zetas)
    zeta = zetas(zi);  % Set current zeta value
    fprintf('zeta = %.1f \n',zeta)

    Or = zeros(length(strength), length(ini_ind));  % Initialize order parameter matrix
    Or_std = zeros(length(strength), length(ini_ind));  % Initialize std deviation matrix
    ACF = zeros(length(strength), length(ini_ind));  % Initialize average autocorrelation matrix
    Or_t = cell(length(ini_ind),1);  % Initialize Order parameter cell
    Finals = cell(length(ini_ind),1); % Initialize final phases cell

    if exist(sprintf('%sOrt_zeta=%.1f.mat',savepath, zeta),'file')
        continue
    end
    tic
    % Loop over different initial conditions
    parfor iter_i = 1:length(ini_ind)
        % Initial frequencies and phases 
        W = W_sets(iter_i,:);
        initial0 = initials(iter_i,:);
        
        % Initialize Matrices
        Or_t_temp= zeros((T(end)-Transient_time)/dt+1,length(strength))  % order parameter time series
        Or_temp = zeros(length(strength), 1);  % average order parameter
        Or_std_temp = zeros(length(strength), 1);  % order parameter std deviation
        ACF_temp = zeros(length(strength), 1);  % average autocorrelation
        finals=zeros(length(strength),size(initial0,2)); % final phases

        % Loop over different coupling strengths
        for si = 1:length(strength)
            if si==1 % if si==1, random initial phases
                initial=initial0; 
            else % if si==1, final phases of the last coupling strength
                initial=finals(si-1,:); 
            end

            % Set the coupling matrix based on current strength
            C = MAT .* strength(si);

            % Simulate the Stuart-Landau model with distance delay and adaptive feedback
            [t,Z] = IE_stuartlandau_distdelay_stim_af(C, dist, speed, noise, W, T, zeta, initial);
            Z=Z(1+Transient_time/dt:end,:); % remove transient time results
            finals(si,:)=Z(end,:);

            % Calculate the global order parameter and its standard deviation
            [or,or_std,or_t] = OrderParameter_Comp(Z);
            Or_temp(si, 1) = or;
            Or_std_temp(si, 1) = or_std;
            Or_t_temp(:,si)=or_t(:);
            acf_tmp = autocorr(or_t,'NumLags',tau/dt); % autocorrelation function
            ACF_temp(si, 1)=acf_tmp(end); % autocorrelation function at tau
        end
        % put temp variables to global variables
        Or(:,iter_i)=Or_temp;
        Or_std(:,iter_i)=Or_std_temp;
        ACF(:,iter_i)=ACF_temp;
        Or_t{iter_i}=Or_t_temp;
        Finals{iter_i}=finals;
    end
    toc
    PCF=size(MAT,1)*Or_std.^2; % Pair Correlation Function

    % save results
    save(sprintf('%sCriticalPoints_zeta=%.1f.mat',savepath, zeta), ...
        'Or', 'Or_std', 'ACF', 'PCF', 'W_sets', 'Finals', 'strength','-v7.3');
    save(sprintf('%sOrt_zeta=%.1f.mat',savepath, zeta),'Or_t','-v7.3');

end

%% KACF, KPCF calculation
fprintf('Calculating KACF, KPCF \n')
KACF_at_critical=zeros(length(zetas),length(ini_ind));
KPCF_at_critical=zeros(length(zetas),length(ini_ind));
for zi = 1:length(zetas)
    zeta = zetas(zi);  % Set current zeta value
    fprintf('zeta = %.1f \n',zeta)

    if exist(sprintf('%sOrt_at_critical_zeta=%.1f.mat',savepath,zeta),'file')

        load(sprintf('%sOrt_at_critical_zeta=%.1f.mat',savepath,zeta),'KACF_temp','KPCF_temp')

    else % calculate KACF and KPCF at critical point

        % Initialize Matrices
        KACF_temp = NaN(length(ini_ind), 1);  % Kurtosis of autocorrelation function
        KPCF_temp = NaN(length(ini_ind), 1);  % Kurtosis of pair correlation function
        Or_t_critical = cell(length(ini_ind),1); % order parameter time series cell
        load(sprintf('%sCriticalPoints_zeta=%.1f.mat',savepath, zeta),'PCF','W_sets','Finals','strength');
        tic
        % Loop over different initial conditions
        parfor iter_i = 1:length(ini_ind)

            % load needed parameters
            finals=Finals{iter_i}
            pcf=PCF(:,iter_i);
            W=W_sets(iter_i,:);

            [~,Max_Ind]=max(pcf); % Critical point is defined as the location of PCF peak
            Max_Ind=max(Max_Ind-shift_peak,1); % Shift 1 coupling strength to avoid discrete dynamics
            C = MAT .* strength(Max_Ind); % coupling strength at critical point

            % if 300sec calculated file is not exist, simulate it
            initial=finals(max(Max_Ind-2,1),:); % initial phase = final phases of previous coupling strength
            [t,Z] = IE_stuartlandau_distdelay_stim_af(C, dist, speed, noise, W, T2, zeta, initial);
            Z=Z(1+Transient_time/dt:end,:); % remove transient time results

            % Calculate the global order parameter and its standard deviation
            [or,or_std,or_t] = OrderParameter_Comp(Z);
            Or_t_critical{iter_i}=or_t;

            % Calculate kurtosis of autocorrelation (tau=0.2sec, windowsize=1sec)
            [KACF_temp(iter_i)] = Ort2KACF(or_t,tau/dt,win_size/dt,overlap);

            % Calculate pair correlation (windowsize=1sec)
            [KPCF_temp(iter_i)] = Ort2KPCF(or_t,size(MAT,1),win_size/dt,overlap);

        end
        save(sprintf('%sOrt_at_critical_zeta=%.1f.mat',savepath,zeta),'Or_t_critical','KACF_temp','KPCF_temp')
        toc
    end
    % put temp variables to global matrix
    KACF_at_critical(zi,:)=KACF_temp;
    KPCF_at_critical(zi,:)=KPCF_temp;
end

% draw figure
figure('Position',[50 50 800 400])
plot_name={'KACF','KPCF'};
for pi=1:length(plot_name)
    subplot(1,2,pi)
    eval(sprintf('plot_data=%s_at_critical;',plot_name{pi}))
    bar(mean(plot_data,2,'omitnan'),'FaceColor',[0.5 0.5 0.5]);
    hold on
    errorbar(1:length(zetas),mean(plot_data,2,'omitnan'),std(plot_data,[],2,'omitnan')/sqrt(sum(~isnan(plot_data),2)),...
        'LineStyle','none','Color',[0.8 0.2 0.2],'LineWidth',0.5)
    ylabel(plot_name{pi})
    xlabel('Z')
    set(gca,'XTickLabel',zetas,'fontsize',12)
end
exportgraphics(gcf,sprintf('%sKACF and KPCF trends.png',savepath),'Resolution',300)
save(sprintf('%sKACF_KPCF_at_critical.mat',savepath),'KACF_at_critical','KPCF_at_critical','zetas')