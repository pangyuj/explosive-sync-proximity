function [S, N2, CCDF] = fss_weighted_network(W, densities)
% Finite-size scaling analysis on a weighted network by percolation
%
% Syntax:
%   [S, N2, CCDF] = fss_weighted_network(W, densities)
%
% Inputs:
%   W          - NxN weighted adjacency matrix (symmetric, can have self-loops but diagonal is ignored)
%   densities  - Array of link densities to scan (e.g., 0.001:0.001:0.06)
%
% Outputs:
%   S     - Structure with fields:
%           .densities - array of link densities scanned
%           .chi       - average susceptibility (s>=2 only, giant excluded) vs density
%           .chi_all   - average susceptibility (all finite, giant excluded) vs density
%   N2    - Structure with fields:
%           .densities - array of link densities scanned
%           .sizes        - average second largest cluster size vs density
%   CCDF  - Structure with fields:
%           .sizes     - unique cluster sizes (excluding giant)
%           .ccdf      - complementary cumulative distribution function values
%           .slope     - power-law exponent fitted over 25-75% quantile range

% Prepare network
N = size(W, 1);
W(1:N+1:end) = 0; % zero out diagonal

% Initialize output arrays
curve_N2 = zeros(numel(densities), 1);
curve_S_all = zeros(numel(densities), 1);
curve_S_ge2 = zeros(numel(densities), 1);

% Pool all finite cluster sizes (excluding giant) for CCDF
sizes_pool = [];

%% Main percolation loop
for di = 1:numel(densities)
    p_dens = densities(di);

    % Build adjacency by selecting EXACTLY top-L edges (tie-safe)
    num_possible = N * (N - 1) / 2;
    L_keep = max(0, min(num_possible, round(p_dens * num_possible)));
    A = adj_by_topL(W, L_keep);

    % Get component sizes
    sizes = getComponentSizes(A);
    sizesExcl = sizes(2:end); % exclude giant

    % N2 (second largest cluster)
    if numel(sizes) >= 2
        curve_N2(di) = sizes(2);
    else
        curve_N2(di) = 0;
    end

    % Susceptibility: chi_all (includes singletons), chi (excludes singletons)
    curve_S_all(di) = susceptibility_from_ns_min(sizesExcl, N, 1);
    curve_S_ge2(di) = susceptibility_from_ns_min(sizesExcl, N, 2);

    % Collect finite cluster sizes for CCDF
    if ~isempty(sizesExcl)
        sizes_pool = [sizes_pool; sizesExcl(:)];
    end

end

%% Compute CCDF
[ccdf_sizes, ccdf_vals] = compute_ccdf(sizes_pool);

% Fit power-law slope over central quantiles (25-75%)
ccdf_slope = NaN;
if ~isempty(ccdf_sizes) && numel(ccdf_sizes) >= 5
    lo = quantile(ccdf_sizes, 0.25);
    hi = quantile(ccdf_sizes, 0.75);
    mid = ccdf_sizes >= lo & ccdf_sizes <= hi & ccdf_vals > 0;

    if nnz(mid) >= 5
        x = log10(ccdf_sizes(mid));
        y = log10(ccdf_vals(mid));
    else
        valid = ccdf_vals > 0;
        x = log10(ccdf_sizes(valid));
        y = log10(ccdf_vals(valid));
    end

    if numel(x) >= 2
        poly_fit = polyfit(x, y, 1);
        ccdf_slope = poly_fit(1);
    end
end

%% Package outputs
S.densities = densities(:);
S.chi = curve_S_ge2(:);
S.chi_all = curve_S_all(:);

N2.densities = densities(:);
N2.sizes = curve_N2(:);

CCDF.sizes = ccdf_sizes(:);
CCDF.ccdf = ccdf_vals(:);
CCDF.slope = ccdf_slope;

end

%% ================== LOCAL HELPER FUNCTIONS ==================

function A = adj_by_topL(W, L_keep)
% Build adjacency by selecting EXACTLY top-L off-diagonal entries (tie-safe)
N = size(W, 1);
maskU = triu(true(N), 1);
vals = W(maskU);

[~, ord] = sort(vals, 'descend');
L_keep = max(0, min(L_keep, numel(ord)));

sel = false(size(vals));
if L_keep > 0
    sel(ord(1:L_keep)) = true;
end

A = zeros(N, N);
A(maskU) = sel;
A = A + A.'; % symmetric, zero diagonal
end

function sizes_sorted = getComponentSizes(A)
% Returns component sizes sorted descending
G = graph(A);
comp = conncomp(G);
K = max(comp);
sizes = accumarray(comp.', 1, [K, 1]);
sizes_sorted = sort(sizes, 'descend');
end

function s = susceptibility_from_ns_min(sizesExcl, N, minSize)
% Susceptibility: giant-excluded, only clusters with size >= minSize
% Formula: s = sum(s^2 * n_s) / sum(s * n_s), where n_s = N_s / N
if nargin < 3, minSize = 1; end
if isempty(sizesExcl), s = 0; return; end

filt = sizesExcl(sizesExcl >= minSize);
if isempty(filt), s = 0; return; end

[uniq, ~, ic] = unique(filt);
counts = accumarray(ic, 1, [numel(uniq), 1]);
ns = counts / N; % n_s = N_s / N

num = sum((uniq.^2) .* ns);
den = sum(uniq .* ns);

if den <= 0
    s = 0;
else
    s = num / den;
end
end

function [sizes_unique, ccdf_vals] = compute_ccdf(sizes_pool)
% Compute complementary cumulative distribution function (CCDF)
if isempty(sizes_pool)
    sizes_unique = [];
    ccdf_vals = [];
    return;
end

sizes_sorted = sort(sizes_pool(:));
sizes_unique = unique(sizes_sorted);

if isempty(sizes_unique)
    ccdf_vals = [];
    return;
end

ccdf_vals = zeros(size(sizes_unique));
for k = 1:numel(sizes_unique)
    ccdf_vals(k) = mean(sizes_sorted >= sizes_unique(k));
end
end
