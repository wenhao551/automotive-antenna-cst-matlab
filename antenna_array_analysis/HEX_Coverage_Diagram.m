clc; clear; close all;

%% =========================================================
%  Build the element radiation pattern G(theta) using HEX_Gain_Phi0
%  First 3 columns in the Excel file:
%   1: Theta [deg]
%   2: Phi   [deg]
%   3: Abs(Dir) [dBi]
%% =========================================================
T = readtable('HEX_Gain_Phi0.xlsx');

theta_all = T{:,1};
phi_all   = T{:,2};
dir_dBi   = T{:,3};

% Fix phi = 0 degrees, and only keep theta = 0~180 degrees
phi0 = 0;
mask = (abs(phi_all - phi0) < 1e-6) & (theta_all >= 0) & (theta_all <= 180);

theta_deg = theta_all(mask);
G_dBi     = dir_dBi(mask);

% dBi -> linear power gain
G_lin = 10.^(G_dBi/10);

% Sort by theta
[theta_deg, idx_sort] = sort(theta_deg);
G_lin = G_lin(idx_sort);

% Remove duplicate theta values to avoid interp1 errors
[theta_deg, ia] = unique(theta_deg, 'stable');
G_lin = G_lin(ia);

% Element pattern interpolation function
% Input angle: geometric angle a_rad relative to the +y axis
% Use |a| in degrees to query the theta cut, assuming left-right symmetry
element_gain = @(a_rad) interp1(theta_deg, G_lin, abs(rad2deg(a_rad)), 'linear', 0);

%% =========================================================
%  Parameters
%% =========================================================
lambda  = 0.012;           % wavelength (m)
k0      = 2*pi/lambda;     % wavenumber
d       = lambda/2;        % element spacing

N       = 20;              % number of elements
Pt_elem = 100;             % Transmit power per element (W)
Pt_tot  = N * Pt_elem;     % Total transmit power (W)
sigma   = 1;               % RCS of each sample point = 1

%% =========================================================
%  Array geometry: along the x-axis and symmetric about the centre
%% =========================================================
idx   = (-(N-1)/2 : (N-1)/2);
ant_x = idx * d;
ant_y = zeros(size(ant_x));

%% =========================================================
%  x-y sample grid
%% =========================================================
x_range = -3:0.05:3;
y_range =  0:0.05:6;
[X, Y] = meshgrid(x_range, y_range);

% Distance from each sample point to the origin, used in the radar equation
R0 = sqrt(X.^2 + Y.^2);
R0(R0 == 0) = eps;

%% =========================================================
%  Scan all steering angles and keep the maximum Pr at each point
%% =========================================================
steer_deg_list = -90:1:90;      % Can be changed to -90:2:90 to speed up

Pr_max     = -inf(size(X));     % Maximum Pr kept for each point
G_best     = nan(size(X));      % Gain when maximum Pr occurs
best_steer = nan(size(X));      % Steering angle when maximum Pr occurs

for s = 1:length(steer_deg_list)

    steer_deg   = steer_deg_list(s);
    theta_steer = deg2rad(steer_deg);

    % Total field under the current steering angle
    Re = zeros(size(X));
    Im = zeros(size(X));

    %% ---------- Element-by-element geometric calculation ----------
    for m = 1:N

        % Relative position from element m to each sample point
        dx = X - ant_x(m);
        dy = Y - ant_y(m);

        % Distance from element m to each sample point
        r_m = sqrt(dx.^2 + dy.^2);
        r_m(r_m == 0) = eps;

        % Direction angle from element m to sample point, relative to +y axis
        a_m = atan2(dx, dy);

        % Look up the element gain in this direction from the HEX table
        Gm = element_gain(a_m);

        % Power density from one element at this point
        Pden_m = Gm .* Pt_elem ./ (4*pi*r_m.^2);

        % Field amplitude
        Em = sqrt(Pden_m);

        % Propagation phase
        phi_prop = k0 * r_m;

        % Feed phase using the actual x-coordinate of the element
        % The steering angle is relative to the +y axis
        phi_feed = -k0 * ant_x(m) * sin(theta_steer);

        % Total phase
        phi_tot = phi_prop + phi_feed;

        % Coherent summation
        Re = Re + Em .* cos(phi_tot);
        Im = Im + Em .* sin(phi_tot);
    end

    % Total power density under the current steering angle
    Pden_arr = hypot(Re, Im).^2;

    % Equivalent gain relative to the total isotropic reference
    Pden_iso = Pt_tot ./ (4*pi*R0.^2);
    G_now = Pden_arr ./ Pden_iso;
    G_now(~isfinite(G_now) | G_now <= 0) = NaN;

    % Radar equation
    Pr_now = Pt_tot .* (G_now.^2) .* lambda.^2 .* sigma ./ ((4*pi)^3 .* R0.^4);
    Pr_now(~isfinite(Pr_now) | Pr_now <= 0) = -inf;

    % Keep the maximum Pr for each sample point
    update_mask = Pr_now > Pr_max;
    Pr_max(update_mask)     = Pr_now(update_mask);
    G_best(update_mask)     = G_now(update_mask);
    best_steer(update_mask) = steer_deg;
end

%% =========================================================
%  Convert to dB
%% =========================================================
Pr_max(Pr_max <= 0 | ~isfinite(Pr_max)) = NaN;
Pr_max_dBW = 10*log10(Pr_max);

G_best(G_best <= 0 | ~isfinite(G_best)) = NaN;
G_best_dB = 10*log10(G_best);

%% =========================================================
%  Display strategy:
%  1) Keep the near-origin region visible
%  2) Exclude a very small region near the origin when setting colour limits
%  3) Clip values above the upper limit to prevent colourbar distortion
%% =========================================================
Pr_plot = Pr_max_dBW;
G_plot  = G_best_dB;
S_plot  = best_steer;

% Ignore the near-origin region only when calculating the colourbar range
r_exclude_for_scale = 0.15;   % Can be changed to 0.10 / 0.20 / 0.30

valid_pr = Pr_max_dBW(isfinite(Pr_max_dBW) & (R0 >= r_exclude_for_scale));
valid_g  = G_best_dB(isfinite(G_best_dB) & (R0 >= r_exclude_for_scale));

% ---------- Manual percentile calculation to avoid using prctile ----------
valid_pr = sort(valid_pr(:));
valid_g  = sort(valid_g(:));

if isempty(valid_pr)
    error('valid_pr is empty. Check the data or r_exclude_for_scale.');
end
if isempty(valid_g)
    error('valid_g is empty. Check the data or r_exclude_for_scale.');
end

idx_pr_hi = max(1, min(numel(valid_pr), round(0.995 * numel(valid_pr))));
idx_g_hi  = max(1, min(numel(valid_g),  round(0.995 * numel(valid_g))));

pr_hi = valid_pr(idx_pr_hi);
pr_lo = pr_hi - 25;      % Can be changed to 30 for a wider range

g_hi = valid_g(idx_g_hi);
g_lo = g_hi - 20;

% Clip values above the display upper limit
Pr_plot(Pr_plot > pr_hi) = pr_hi;
G_plot(G_plot > g_hi) = g_hi;

%% =========================================================
%  Figure 1: maximum Pr at each point
%% =========================================================
figure;
imagesc(x_range, y_range, Pr_plot);
axis xy; axis equal; box on;
colormap(jet);
cb = colorbar;
ylabel(cb, 'Max P_r (dBW)');
caxis([pr_lo pr_hi]);

title('Coverage Diagram');
xlabel('Lateral position, x (m)');
ylabel('Longitudinal position, y (m)');