clc; clear; close all;

%% =========================================================
%  Build the element radiation pattern G(theta)
%  using CIR_Gain_Phi90 data exported from CST.
%
%  Excel file CIR_Gain_Phi90.xlsx columns:
%   1: Theta [deg]
%   2: Phi   [deg]
%   3: Abs(Gain) [dBi]
%
%  Main idea:
%   1) Extract the phi = 90 degree cut from the CST radiation pattern.
%   2) Convert gain from dBi to linear scale.
%   3) Build an interpolation function element_gain(a_rad).
%   4) Use the real element pattern in array field summation and Pr calculation.
%% =========================================================
T = readtable('CIR_Gain_Phi90.xlsx');

% Read theta, phi, and gain in dBi from the Excel file
theta_all = T{:,1};
phi_all   = T{:,2};
dir_dBi   = T{:,3};

% Select the phi = 90 degree cut
phi0 = 90;
mask = (abs(phi_all - phi0) < 1e-6);
% Use mask to select data where phi = 90 degrees

theta_deg = theta_all(mask);
G_dBi     = dir_dBi(mask);

% Convert dBi gain to linear power gain
G_lin = 10.^(G_dBi/10);

% Sort theta values for interpolation
[theta_deg, idx_sort] = sort(theta_deg);% Sort theta in ascending order and record the original indices
G_lin = G_lin(idx_sort);% Reorder gain values using the same index order

% ==========================================================
% Element pattern interpolation function
%
% Input:
%   a_rad : observation angle relative to the +y axis, in radians
%
% Output:
%   linear element gain in that direction
%
% abs(rad2deg(a_rad)) is used because the pattern is treated
% as symmetric about the main axis.
%
% interp1(..., 'linear', 0):
%   - uses linear interpolation between known data points
%   - returns 0 if the query angle is outside the data range
% ==========================================================
% Function input: direction angle a_rad in radians. Output: element gain
element_gain = @(a_rad) interp1(theta_deg, G_lin, abs(rad2deg(a_rad)), 'linear', 0);

%% =========================================================
%  Basic parameter settings
%% =========================================================

% Operating wavelength (m)
lambda  = 0.012;

% Wavenumber: phase accumulated per metre of propagation
k0      = 2*pi/lambda;

% Element spacing, set to lambda/2
d       = lambda/2;

% Number of array elements
N       = 20;

% Transmit power of each element (W)
Pt_elem = 100;

% Total array transmit power = number of elements x element power
Pt_tot  = N * Pt_elem;

% Radar Cross Section (RCS) of each sample point
% Here sigma = 1 is used for relative comparison
sigma   = 1;

%% =========================================================
%  Array geometry: 20 elements placed along the x-axis,
%  symmetric about the origin
idx   = (-(N-1)/2 : (N-1)/2);  % Centred element index
ant_x = idx * d;               % x-coordinate of each element
ant_y = zeros(size(ant_x));    % All elements are located at y = 0

%% =========================================================
%  Build spatial sample grid
%
%  x: -10 m to 10 m
%  y:  0 m to 10 m
%
%  This represents the front half-plane of the array.
%% =========================================================
x_range = -3:0.05:3;
y_range =   0:0.05:6;
[X, Y] = meshgrid(x_range, y_range);

% Distance from each sample point to the origin
R0 = sqrt(X.^2 + Y.^2);

% Avoid division by zero at the origin
R0(R0 == 0) = eps;%

%% =========================================================
%  Scan all steering angles and keep the maximum Pr
%  at each sample point.
%
%  For each steering angle:
%   - calculate the full field map;
%   - compare the current Pr with the previous maximum;
%   - update the maximum Pr, corresponding G, and steering angle if larger.
%% =========================================================
steer_deg_list = -90:1:90;      % Steering angles relative to +y axis. +90 degrees = right, -90 degrees = left

% Initial values:
% Pr_max     : maximum received power at each point
% G_best     : equivalent array gain when maximum Pr occurs
% best_steer : steering angle when maximum Pr occurs
Pr_max     = -inf(size(X));
G_best     = nan(size(X));
best_steer = nan(size(X));

for s = 1:length(steer_deg_list)% Outer loop: calculate the full map for each steering angle

    % Current steering angle in degrees
    steer_deg   = steer_deg_list(s);

    % Convert to radians for phase calculation
    theta_steer = deg2rad(steer_deg);

    % Initialize real and imaginary parts of the total array field
    Re = zeros(size(X));
    Im = zeros(size(X));

    %% ---------- Element-by-element coherent field summation ----------
    for m = 1:N % Inner loop: calculate the contribution from each element

        %% Distance from element m to each sample point
        % Relative coordinates from element m to each sample point
        dx = X - ant_x(m);% ant_x(m): x-coordinate of element m
        dy = Y - ant_y(m);% ant_y(m): y-coordinate of element m

        % Propagation distance from element m to each sample point
        r_m = sqrt(dx.^2 + dy.^2);
        r_m(r_m == 0) = eps;   % Avoid division by zero

        % Direction angle from element m to each sample point
        % The angle is defined relative to the +y axis
        a_m = atan2(dx, dy);

        % Look up the element gain from the real element pattern
        Gm = element_gain(a_m);

        % --------------------------------------------------
        % Power density generated by one element at the sample point:
        %
        %   Pden = G * Pt / (4*pi*r^2)
        %
        % The element transmit power is weighted by the element gain
        % and attenuated by spherical spreading.
        % --------------------------------------------------
        Pden_m = Gm .* Pt_elem ./ (4*pi*r_m.^2);

        % Convert power density to field amplitude
        Em = sqrt(Pden_m);

        % Propagation phase from element m to each sample point
        phi_prop = k0 * r_m; % Propagation phase determined by distance

        % Feed phase for beam steering
        % This phase makes the main lobe point toward the steering angle
        phi_feed = -k0 * ant_x(m) * sin(theta_steer);% Feed phase set by the selected steering angle

        % Total phase = propagation phase + feed phase
        phi_tot = phi_prop + phi_feed;

        % --------------------------------------------------
        % Coherent summation:
        % Treat each element field as a complex field E * exp(j*phi)
        % and add real and imaginary parts separately.
        % --------------------------------------------------
        Re = Re + Em .* cos(phi_tot);
        Im = Im + Em .* sin(phi_tot);
    end

    % ======================================================
    % Total power density for the current steering angle
   
    Pden_arr = hypot(Re, Im).^2;% Equivalent to Re.^2 + Im.^2;

    % ======================================================
    % Isotropic reference power density with the same total power.
    %
    % If the total array power Pt_tot were radiated by an ideal
    % isotropic radiator, the power density at distance R0 would be:
    %   Pden_iso = Pt_tot / (4*pi*R0^2)
    %
    % The ratio between the actual array power density and this
    % isotropic reference gives the equivalent array gain G_now.
    % ======================================================
    Pden_iso = Pt_tot ./ (4*pi*R0.^2);% Pt_tot: total array transmit power
    G_now = Pden_arr ./ Pden_iso;

    % Remove invalid or non-positive values
    G_now(~isfinite(G_now) | G_now <= 0) = NaN;

    % ======================================================
    % Calculate received power Pr using the radar equation:
    %
    %   Pr = Pt * G^2 * lambda^2 * sigma / ((4*pi)^3 * R^4)
    %
    % G_now is the equivalent array gain at each sample point.
    % R0 is the distance from the array reference centre.
    % sigma is the assumed RCS.
    %
    % This is used as a relative spatial Pr estimate for
    % coverage diagram visualisation and antenna comparison.
    % ======================================================
    % Calculate received power Pr using the radar equation
    Pr_now = Pt_tot .* (G_now.^2) .* lambda.^2 .* sigma ./ ((4*pi)^3 .* R0.^4);

    % Set invalid or non-positive values to -Inf for max comparison
    Pr_now(~isfinite(Pr_now) | Pr_now <= 0) = -inf;

    % ======================================================
    % Keep only the maximum Pr at each sample point.
    % Also record the corresponding G and steering angle.
    % ======================================================
    update_mask = Pr_now > Pr_max;
    Pr_max(update_mask)     = Pr_now(update_mask);
    G_best(update_mask)     = G_now(update_mask);
    best_steer(update_mask) = steer_deg;
end

%% =========================================================
%  Convert results to dB scale for plotting
%% =========================================================

% Convert Pr to dBW
Pr_max(Pr_max <= 0 | ~isfinite(Pr_max)) = NaN;
Pr_max_dBW = 10*log10(Pr_max);

% Convert G to dB
G_best(G_best <= 0 | ~isfinite(G_best)) = NaN;
G_best_dB = 10*log10(G_best);

%% =========================================================
%  Display strategy
%
%  Values near the origin can become very large because the
%  distance is very small. This can make the colourbar range too wide.
%
%  Therefore:
%   1) keep the near-origin region in the plot;
%   2) exclude a small region near the origin when estimating colour limits;
%   3) clip values above the upper colour limit.
%% =========================================================
Pr_plot = Pr_max_dBW;
G_plot  = G_best_dB;
S_plot  = best_steer;

% Exclude points too close to the origin when calculating colour limits
r_exclude_for_scale = 0.15;   % Can be changed to 0.10 / 0.20 / 0.30 if needed

% Valid values used to calculate the Pr colour range
valid_pr = Pr_max_dBW(isfinite(Pr_max_dBW) & (R0 >= r_exclude_for_scale));

% Valid values used to calculate the G colour range
valid_g  = G_best_dB(isfinite(G_best_dB) & (R0 >= r_exclude_for_scale));

% ==========================================================
% Manual percentile calculation.
%
% prctile is not used to avoid dependency issues in some
% MATLAB environments.
% ==========================================================
valid_pr = sort(valid_pr(:));% Sort all valid values in ascending order for percentile calculation
valid_g  = sort(valid_g(:));% Sort all valid values in ascending order for percentile calculation

% Use the 99.5th percentile as the upper display limit
idx_pr_hi = max(1, min(numel(valid_pr), round(0.995 * numel(valid_pr))));
idx_g_hi  = max(1, min(numel(valid_g),  round(0.995 * numel(valid_g))));

pr_hi = valid_pr(idx_pr_hi);
g_hi  = valid_g(idx_g_hi);

% Lower limits are defined by fixed dynamic ranges
pr_lo = pr_hi - 25;      % 25 dB dynamic range for Pr plot
g_lo  = g_hi  - 20;      % 20 dB dynamic range for G plot

% Clip values above the display upper limit
Pr_plot(Pr_plot > pr_hi) = pr_hi;
G_plot(G_plot > g_hi) = g_hi;

%% =========================================================
%  Figure 1: maximum received power Pr at each sample point
%
%  x-axis: x (m)
%  y-axis: y (m)
%  colour: maximum Pr over all steering angles, in dBW
%% =========================================================
figure;
imagesc(x_range, y_range, Pr_plot);
axis xy;        % Keep the y-axis direction upward
axis equal;     % Keep equal x/y scale
box on;

colormap(jet);
cb = colorbar;
ylabel(cb, 'Max P_r (dBW)');

% Set colourbar range
caxis([pr_lo pr_hi]);

title('Coverage Diagram');
xlabel('Lateral position, x (m)');
ylabel('Longitudinal position, y (m)');