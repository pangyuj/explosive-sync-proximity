% Example usage of fss_weighted_network

% Create a sample weighted network
N = 1000;
W = rand(N);
W = (W + W') / 2; % make symmetric
W(1:N+1:end) = 0; % zero out diagonal

% Define density range
densities = 0.0001:0.0001:0.006;

% Perform FSS analysis
[S, N2, CCDF] = fss_weighted_network(W, densities);

% Create figure with subplots
figure('Position', [100 100 1200 400]);

% Subplot 1: Susceptibility (chi)
subplot(1, 3, 1);
plot(S.densities, S.chi, 'k-', 'LineWidth', 1.5);
hold on;
[chi_max, idx_max] = max(S.chi);
dens_max = S.densities(idx_max);
plot(dens_max, chi_max, 'r.', 'MarkerSize', 15, 'LineWidth', 2);
xlabel('Link Density');
ylabel('Susceptibility \chi (s\geq2)');
grid on;
title('Susceptibility vs Density');
% Subplot 2: Second Largest Cluster
subplot(1, 3, 2);
plot(N2.densities, N2.sizes, 'k-', 'LineWidth', 1.5);
hold on
plot(dens_max, max(N2.sizes), 'r.', 'MarkerSize', 15, 'LineWidth', 2);
xlabel('Link Density');
ylabel('N_2');
grid on;
title('Second Largest Cluster vs Density');

% Subplot 3: CCDF (log-log)
subplot(1, 3, 3);
loglog(CCDF.sizes, CCDF.ccdf, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
hold on;
if ~isnan(CCDF.slope)
    lo = quantile(CCDF.sizes, 0.25);
    hi = quantile(CCDF.sizes, 0.75);
    mid = CCDF.sizes >= lo & CCDF.sizes <= hi;
    
    x_fit = [lo, hi];
    y_fit = 10.^(log10(CCDF.ccdf(find(mid, 1))) + CCDF.slope * (log10(x_fit) - log10(lo)));
    loglog(x_fit, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('slope: %.2f', CCDF.slope));
    legend;
end
xlabel('Cluster Size');
ylabel('CCDF');
grid on;
title('Cluster Size Distribution');
