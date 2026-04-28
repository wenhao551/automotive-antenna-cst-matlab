clc; clear; close all;

%% ---------- Build element radiation pattern G(phi) using Gain_Theta25 ----------
% First 3 columns in the Excel file:
%   1: Theta [deg]
%   2: Phi   [deg]
%   3: Abs(Dir) [dBi]

T = readtable('SQU_Gain_Theta25.xlsx');   % If there is no extension, change to 'Gain_Theta25'

theta_all = T{:,1};
phi_all   = T{:,2};
dir_dBi   = T{:,3};

% Fix theta = 25 degrees, and only keep phi = 0~180 degrees
theta0 = 25;      % Change this value if another cut is needed
mask = (abs(theta_all - theta0) < 1e-6) & (phi_all >= 0) & (phi_all <= 180);

phi_deg = phi_all(mask);    % 0~180 degrees
G_dBi   = dir_dBi(mask);

% dBi -> linear gain
G_lin = 10.^(G_dBi/10);

% Sort by phi
[phi_deg, idx] = sort(phi_deg);
G_lin          = G_lin(idx);

% Element pattern:
% a_rad is the geometric angle from the element to the observation point,
% measured relative to the +y axis, in radians.
% Use |angle| -> [0,180] for table lookup, so the pattern is treated as left-right symmetric.
element_gain = @(a_rad) ...
    interp1(phi_deg, G_lin, ...
            abs(rad2deg(a_rad)), ...   % Absolute value
            'linear', 0);              % Return 0 outside the data range

%% ---------- Basic parameters ----------
lambda = 0.012;          % Wavelength (m)
d = lambda/2;            % Element spacing (m)
N = 20;                  % Number of elements
Pt = 100;                % Transmit power per element (W)

%% ---------- Array geometry: along the x-axis and symmetric about the origin ----------
idx = (-(N-1)/2 : (N-1)/2);   % -9.5:9.5 (length N)
ant_x = idx * d;              % Element x-coordinates
ant_y = zeros(size(ant_x));   % All y = 0

%% ---------- Observation grid ----------
x_range = -10:0.05:10;        % x range (m)
y_range =   0:0.05:20;        % y range (m)
[X, Y] = meshgrid(x_range, y_range);

%% ---------- Reference: isotropic single-antenna power density ----------
dx1 = X;
dy1 = Y;
r1  = sqrt(dx1.^2 + dy1.^2);
r1(r1 == 0) = eps;            % Avoid division by zero

Pden_iso = Pt ./ (4*pi*r1.^2);    % Isotropic reference, without element pattern

%% ---------- 20-element array power density using complex-amplitude summation ----------
Re = zeros(size(X));
Im = zeros(size(X));

for k = 1:N
    dx = X - ant_x(k);
    dy = Y - ant_y(k);
    r  = sqrt(dx.^2 + dy.^2);
    r(r == 0) = eps;

    % Geometric radiation angle a relative to +y, approximately within [-90, 90] degrees
    a = atan2(dx, dy);

    % Element pattern power gain from Gain_Theta25 interpolation
    G_elem = element_gain(a);

    % Power density from one element in this direction: Pd_k = G_elem * Pt / (4*pi*r^2)
    Pden_k = G_elem .* Pt ./ (4*pi*r.^2);

    % Field amplitude is proportional to sqrt(power density)
    V_k = sqrt(Pden_k);

    % Propagation phase: phi_k = 2*pi*r/lambda
    phi_k = 2*pi*r/lambda;

    % Coherent summation
    Re = Re + V_k .* cos(phi_k);
    Im = Im + V_k .* sin(phi_k);
end

Vabs = hypot(Re, Im);            % Magnitude of total field
Pden_array20 = Vabs.^2;          % Array power density

%% ---------- Total array gain in linear scale and dB, relative to N isotropic elements ----------
G = Pden_array20 ./ (N * Pden_iso);
G(~isfinite(G)) = NaN;
GdB = 10*log10(G);

%% ---------- Visualisation: 2D gain distribution in dB ----------
figure;
imagesc(x_range, y_range, GdB);
axis xy; axis equal; box on;
colormap(jet); colorbar;
clim([-20 25]);               % Adjust if needed
title('Total Gain Map of a 20-Element Linear Array');
xlabel('Lateral position, x (m)');
ylabel('Longitudinal position, y (m)');
cb = colorbar;
ylabel(cb, 'Array gain, G (dB)');