% This script simulates a Stuart-Landau network with distance delay and adaptive feedback.
% Calculates the time to lose and regain criticality.
clear; clc;

% set save path and generate folder if not exist
savepath='results_400iter/';

% simulation parameters
load('gong78');  % Load connectivity (MAT) and distance matrix (Dmat)
dist = Dmat/1000;  % Convert distance to meters (m)
noise = 1;        % Set noise level
speed = 7;        % Speed of signal propagation (m/s)
dt = 0.001;       % Time step for the simulation
T = 0:dt:25;      % Time vector (simulation duration 25 seconds)
Stimulation_time=10; % Stimulation at 10 seconds
Stim_str=40; % Stimulation strength
Stim_dur=5; % Stimulation duration

% Define zeta values (adaptive feedback strength)
zetas = 0:0.5:3;

% Initialize lose & regain time matrices
Time_to_lose_criticality=[];
Time_to_regain_criticality=[];

fprintf('Stim_str = %.1f Stim_dur=%.3f \n',Stim_str,Stim_dur)

% Loop over each zeta value (adaptive feedback strength)
for zi = 1:length(zetas)
    zeta = zetas(zi);  % Set current zeta value
    fprintf('zeta = %.1f \n',zeta)

    % load critical point information for a zeta
    load(sprintf('%sCriticalPoints_zeta=%.1f.mat', savepath, zeta), ...
        'Or', 'Or_std', 'ACF', 'W_sets', 'Finals', 'strength');

    % Generate stimulation time series (u)
    Stim_start=find(Stimulation_time==T);
    Stim_end=Stim_start-1+Stim_dur/dt;
    u=T*0; u(Stim_start:Stim_end)=Stim_str*(1+1i);
    tic
    parfor iter_i = 1:size(W_sets,1)
        % Find critical point with PCF peak
        [~,Max_Ind]=max(Or_std(:,iter_i));
        Max_Ind=Max_Ind-1; % Shift 1 coupling strength to avoid discrete dynamics

        % Get initial conditions
        W = W_sets(iter_i,:); % node frequencies
        finals=Finals{iter_i}; % load final phases
        initial=finals(max(Max_Ind-1,1),:); % initial phase = final phases of previous coupling strength

        % Set the coupling matrix based on coupling strength at critical point
        C = MAT .* strength(Max_Ind);

        % Simulate the Stuart-Landau model with distance delay and adaptive feedback
        [t,Z] = IE_stuartlandau_distdelay_stim_af(C, dist, speed, noise, W, T, zeta, initial, u);

        % Calculate the global order parameter
        [~,~,or_t] = OrderParameter_Comp(Z);

        % Threshold for criticality (3 standard deviations)
        criticality_thr=Or(Max_Ind,iter_i)+3*Or_std(Max_Ind,iter_i);

        % Calculate time to lose criticality
        lose_ti=find(or_t(Stim_start:Stim_end)>criticality_thr); % Time index for losing criticality
        if ~isempty(lose_ti)
            Time_to_lose_criticality(zi,iter_i)=lose_ti(1)*dt; % First time index exceeding threshold
        else
            Time_to_lose_criticality(zi,iter_i)=NaN;
        end

        % Calculate time to regain criticality
        regain_ti=find(or_t(Stim_end+1:end)<criticality_thr); % Time index for regaining criticality
        if ~isempty(regain_ti)
            Time_to_regain_criticality(zi,iter_i)=regain_ti(1)*dt; % First time index below threshold
        else
            Time_to_regain_criticality(zi,iter_i)=NaN;
        end
    end
    toc

end
% save results
save(sprintf('%sTimeToLoseRegainCriticality.mat', savepath), ...
    'Time_to_lose_criticality', 'Time_to_regain_criticality', 'Stim_str', 'Stim_dur', '-v7.3');
%% plot results
load(sprintf('%sTimeToLoseRegainCriticality.mat', savepath))
load(sprintf('%sKACF_KPCF_at_critical.mat',savepath))

figure('position',[50 50 800 400])
% KACF values
x=mean(KACF_at_critical,2,'omitnan');
x_low=std(KACF_at_critical,[],2,'omitnan')./sqrt(sum(~isnan(KACF_at_critical),2));
x_high=std(KACF_at_critical,[],2,'omitnan')./sqrt(sum(~isnan(KACF_at_critical),2));

% Time to lose criticality
subplot(1,2,1)
y=mean(Time_to_lose_criticality,2,'omitnan');
y_low=std(Time_to_lose_criticality,[],2,'omitnan')./sqrt(sum(~isnan(Time_to_lose_criticality),2));
y_high=std(Time_to_lose_criticality,[],2,'omitnan')./sqrt(sum(~isnan(Time_to_lose_criticality),2));
e1=errorbar(x,y,y_low,y_high,x_low,x_high,'LineStyle','none','Color',[0.75 0.1 0],'linewidth',1);
e1.CapSize=0;
xlabel('Kurtosis of autocorrelation function')
ylabel('Time to lose criticality (sec)','Color',[0.75 0.1 0])

% Time to regain criticality
subplot(1,2,2)
y=mean(Time_to_regain_criticality,2,'omitnan');
y_low=std(Time_to_regain_criticality,[],2,'omitnan')./sqrt(sum(~isnan(Time_to_regain_criticality),2));
y_high=std(Time_to_regain_criticality,[],2,'omitnan')./sqrt(sum(~isnan(Time_to_regain_criticality),2));
e1=errorbar(x,y,y_low,y_high,x_low,x_high,'LineStyle','none','Color',[0.0 0.45 0.75],'linewidth',1);
e1.CapSize=0;
xlabel('Kurtosis of autocorrelation function')
ylabel('Time to regain criticality (sec)','Color',[0.0 0.45 0.75])
exportgraphics(gcf,sprintf('%sTime to Lose, Regain Criticality.png',savepath),'Resolution',300)

