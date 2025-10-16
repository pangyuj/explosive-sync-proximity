% This script simulates a Stuart-Landau network with distance delay and adaptive feedback.
% Calculates the time to lose and regain criticality.
% To test robustness on stimulation parameters
clear; clc;

% set save path and generate folder if not exist
savepath='results_400iter/';

% simulation parameters
load('gong78');  % Load connectivity (MAT) and distance matrix (Dmat)
dist = Dmat/1000;  % Convert distance to meters (m)
noise = 1;        % Set noise level
speed = 7;        % Speed of signal propagation (m/s)
dt = 0.001;       % Time step for the simulation
T = 0:dt:40;      % Time vector (simulation duration 40 seconds)
Stimulation_time=10; % Stimulation at 10 seconds
Stimulation_strengths=20:20:100; % Range of stimulation strengths
Stimulation_durations=[0.1 0.2 0.5 1 2 5 10 20]; % Range of stimulation durations

% Define zeta values (adaptive feedback strength)
zetas = 0:3;

% Loop over each zeta value (adaptive feedback strength)
for zi = 1:length(zetas)
    zeta = zetas(zi);  % Set current zeta value
    fprintf('zeta = %.1f \n',zeta)

    % if result file exist, continue to other zeta
    if exist(sprintf('%sTimeToLoseRegainCriticality_PCFpeak_zeta=%.1f.mat', savepath, zeta),'file')
        continue
    end
    
    % load critical point information for a zeta
    load(sprintf('%sCriticalPoints_zeta=%.1f.mat', savepath, zeta), ...
        'Or', 'Or_std', 'ACF', 'W_sets', 'Finals', 'strength');

    % Initialize lose & regain time matrices
    Time_to_lose_criticality=zeros(length(Stimulation_strengths),length(Stimulation_durations),size(W_sets,1));
    Time_to_regain_criticality=zeros(length(Stimulation_strengths),length(Stimulation_durations),size(W_sets,1));
    tic
    for str_i=1:length(Stimulation_strengths)
        Stim_str=Stimulation_strengths(str_i);
        for dur_i=1:length(Stimulation_durations)
            Stim_dur=Stimulation_durations(dur_i);
            fprintf('Stim_str = %.1f Stim_dur=%.3f, ',Stim_str,Stim_dur)
            
            % Generate stimulation time series (u)
            Stim_start=find(Stimulation_time==T);
            Stim_end=Stim_start-1+Stim_dur/dt;
            u=T*0; u(Stim_start:Stim_end)=Stim_str*(1+1i);
            
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
                    Time_to_lose_criticality(str_i,dur_i,iter_i)=lose_ti(1)*dt; % First time index exceeding threshold
                else
                    Time_to_lose_criticality(str_i,dur_i,iter_i)=NaN;
                end

                % Calculate time to regain criticality
                regain_ti=find(or_t(Stim_end+1:end)<criticality_thr); % Time index for regaining criticality
                if ~isempty(regain_ti)
                    Time_to_regain_criticality(str_i,dur_i,iter_i)=regain_ti(1)*dt; % First time index below threshold
                else
                    Time_to_regain_criticality(str_i,dur_i,iter_i)=NaN;
                end
            end
            toc
        end
    end
    % save results
    save(sprintf('%sTimeToLoseRegainCriticality_RobustnessTest_zeta=%.1f.mat', savepath, zeta), ...
        'Time_to_lose_criticality', 'Time_to_regain_criticality', 'Stimulation_strengths', 'Stimulation_durations', '-v7.3');
end

%% plot results
close all;
% Preallocate matrices to store mean and standard error of the mean (SEM) for time to lose and regain criticality
Mean_Time_to_lose = zeros(length(Stimulation_strengths), length(Stimulation_durations), length(zetas));
SEM_Time_to_lose = zeros(length(Stimulation_strengths), length(Stimulation_durations), length(zetas));
Mean_Time_to_regain = zeros(length(Stimulation_strengths), length(Stimulation_durations), length(zetas));
SEM_Time_to_regain = zeros(length(Stimulation_strengths), length(Stimulation_durations), length(zetas));

% Loop over zeta values
for zi = 1:length(zetas)
    zeta = zetas(zi);
    
    % Load the time to lose and regain criticality for the current zeta value
    load(sprintf('%sTimeToLoseRegainCriticality_RobustnessTest_zeta=%.1f.mat',savepath, zeta), ...
        'Time_to_lose_criticality', 'Time_to_regain_criticality', 'Stimulation_strengths', 'Stimulation_durations');
    Time_to_regain_criticality(Time_to_regain_criticality<=dt)=NaN;
    % Compute the mean time to lose and regain criticality
    Mean_Time_to_lose(:, :, zi) = mean(Time_to_lose_criticality, 3, 'omitnan');
    Mean_Time_to_regain(:, :, zi) = mean(Time_to_regain_criticality, 3, 'omitnan');
    
    % Loop over stimulation strengths and durations
    for str_i = 1:length(Stimulation_strengths)
        for dur_i = 1:length(Stimulation_durations)
            % Compute the SEM for time to lose and regain criticality
            SEM_Time_to_lose(str_i, dur_i, zi) = std(Time_to_lose_criticality(str_i, dur_i, :), [], 3, 'omitnan') / ...
                sqrt(sum(~isnan(Time_to_lose_criticality(str_i, dur_i, :)), 3));
            SEM_Time_to_regain(str_i, dur_i, zi) = std(Time_to_regain_criticality(str_i, dur_i, :), [], 3, 'omitnan') / ...
                sqrt(sum(~isnan(Time_to_regain_criticality(str_i, dur_i, :)), 3));
        end
    end
end

% Create a figure to plot the mean and SEM for time to lose criticality
figure('position', [50 50 1200 500]); plot_i = 1;
plot_colors=jet(length(Stimulation_strengths))*0.8;
% Loop over stimulation strengths and durations for plotting
for dur_i = 1:length(Stimulation_durations)
    subplot(2, 4, dur_i);
    for str_i = 1:length(Stimulation_strengths)
          % Plot the mean and SEM for each combination of strength and duration
        plot_i = plot_i + 1;
        plot_mean = 1000*Mean_Time_to_lose(str_i, dur_i, :);
        plot_mse = 1000*SEM_Time_to_lose(str_i, dur_i, :);
        
        % Plot error bars showing mean and SEM
        errorbar(zetas+0.025*(str_i-length(Stimulation_strengths)/2), plot_mean(:), plot_mse(:), 'LineWidth', 1, 'Color',plot_colors(str_i,:));
        hold on
    end
    set(gca,'xticklabel',append('Z=',string(num2str(zetas'))),'yscale','log')
    ylabel('time (ms)')
    grid on
    title(sprintf('duration = %.1f sec', Stimulation_durations(dur_i)));
end
annotation(gcf,"textbox",'LineStyle','none','Position',[0.28,0.965,0.5,0.05], 'HorizontalAlignment','center', 'Fontsize', 14,...
    'String','Time to Lose Criticality')
legend(append('p=',string(num2str(Stimulation_strengths'))),'position',[0.92 0.36 0.069 0.31],'Box','off')
exportgraphics(gcf,sprintf('%sTime to Lose Criticality Robustness test.png',savepath),'Resolution',300)


% Create a figure to plot the mean and SEM for time to regain criticality
figure('position', [50 50 1200 500]); plot_i = 1;
% Loop over stimulation strengths and durations for plotting
for dur_i = 1:length(Stimulation_durations)
    subplot(2, 4, dur_i);
    for str_i = 1:length(Stimulation_strengths)
        % Plot the mean and SEM for each combination of strength and duration
        plot_i = plot_i + 1;
        plot_mean = 1000*Mean_Time_to_regain(str_i, dur_i, :);
        plot_mse = 1000*SEM_Time_to_regain(str_i, dur_i, :);
        
        % Plot error bars showing mean and SEM
        errorbar(zetas+0.025*(str_i-length(Stimulation_strengths)/2), plot_mean(:), plot_mse(:), 'LineWidth', 1, 'Color',plot_colors(str_i,:));
        hold on
    end
    set(gca,'xticklabel',append('Z=',string(num2str(zetas'))),'yscale','log')
    ylabel('time (ms)')
    grid on
    title(sprintf('duration = %.1f sec', Stimulation_durations(dur_i)));
end
annotation(gcf,"textbox",'LineStyle','none','Position',[0.28,0.965,0.5,0.05], 'HorizontalAlignment','center', 'Fontsize', 14,...
    'String','Time to Regain Criticality')
legend(append('p=',string(num2str(Stimulation_strengths'))),'position',[0.92 0.36 0.069 0.31],'Box','off')
exportgraphics(gcf,sprintf('%sTime to Regain Criticality Robustness test.png',savepath),'Resolution',300)