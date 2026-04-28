clc; clear; close all;

%% ---------- Build element radiation pattern G(theta) using CIR_Gain_Phi90 ----------
% First 3 columns in the Excel file:
%   1: Theta [deg]
%   2: Phi   [deg]
%   3: Abs(Gain)[dBi]

T = readtable('CIR_Gain_Phi90.xlsx');

theta_all = T{:,1};
phi_all   = T{:,2};
dir_dBi   = T{:,3};

% Fix phi = 90 degrees, and only keep theta = 0~180 degrees
phi0 = 90;
mask = (abs(phi_all - phi0) < 1e-6) & (theta_all >= 0) & (theta_all <= 180);

theta_deg = theta_all(mask);     % 0~180
G_dBi     = dir_dBi(mask);

% dBi -> linear power gain
G_lin = 10.^(G_dBi/10);

% Sort by theta
[theta_deg, idx] = sort(theta_deg);
G_lin             = G_lin(idx);

% Element pattern interpolation function
% Map the geometric angle a to theta for table lookup
% a_rad: geometric angle relative to +y, in [-pi, pi]
% theta_query = |a| (deg) -> [0,180], assuming left-right symmetry
element_gain = @(a_rad) ...
    interp1(theta_deg, G_lin, ...
            abs(rad2deg(a_rad)), ...
            'linear', 0);

%% ---------- Basic parameters ----------
lambda = 0.012;     % wavelength (m)
d = lambda/2;       % element spacing
N = 20;             % number of elements
Pt = 100;           % power per element (W)

%% ---------- Array geometry: along the x-axis and symmetric about the origin ----------
idx = (-(N-1)/2 : (N-1)/2);
ant_x = idx * d;
ant_y = zeros(size(ant_x));

%% ---------- Observation grid ----------
x_range = -10:0.05:10;
y_range =   0:0.05:20;
[X, Y] = meshgrid(x_range, y_range);

%% ---------- Reference: isotropic single-antenna power density ----------
r1 = hypot(X, Y);
r1(r1==0) = eps;
Pden_iso = Pt ./ (4*pi*r1.^2);

%% ---------- 20-element array power density using complex-amplitude vector summation ----------
Re = zeros(size(X));
Im = zeros(size(X));

for k = 1:N
    dx = X - ant_x(k);
    dy = Y - ant_y(k);

    r = hypot(dx, dy);
    r(r==0) = eps;

    % Geometric angle a relative to +y
    a = atan2(dx, dy);

    % Element pattern from the phi = 90 degree cut
    G_elem = element_gain(a);

    % Single-element power density
    Pden_k = G_elem .* Pt ./ (4*pi*r.^2);

    % Field amplitude is proportional to sqrt(power density)
    V_k = sqrt(Pden_k);

    % Propagation phase
    phi_k = 2*pi*r/lambda;

    % Coherent summation
    Re = Re + V_k .* cos(phi_k);
    Im = Im + V_k .* sin(phi_k);
end

Vabs = hypot(Re, Im);
Pden_array20 = Vabs.^2;

%% ---------- Total array gain relative to N isotropic elements ----------
G = Pden_array20 ./ (N * Pden_iso);
G(~isfinite(G)) = NaN;
GdB = 10*log10(G);

%% ---------- Visualisation ----------
figure;
imagesc(x_range, y_range, GdB);
axis xy; axis equal; box on;
colormap(jet); colorbar;
clim([-20 25]);

title('Total Gain Map of a 20-Element Linear Array');
xlabel('Lateral position, x (m)');
ylabel('Longitudinal position, y (m)');
cb = colorbar;
ylabel(cb, 'Array gain, G (dB)');