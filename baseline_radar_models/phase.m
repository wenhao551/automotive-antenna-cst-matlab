clc; clear; close all;

% ---------- Parameters ----------
lambda = 0.012;          % Wavelength (m)
d = lambda/2;            % Spacing between adjacent antennas
N = 20;                  % Number of antennas
Pt = 1;                  % Transmit power of each antenna (W)

delta_deg = 60;                  % Phase difference between adjacent array elements (degrees)
delta_rad = deg2rad(delta_deg);  % Convert to radians

% Antenna coordinates: symmetric about the y-axis, with no element at the origin when the number is even
idx = (-(N-1)/2 : (N-1)/2);  % -9.5:9.5
ant_x = idx * d;             % x positions (increasing from left to right)
ant_y = zeros(size(ant_x));  % All y=0

% Calculate grid
x_range = -10:0.05:10;       % x-axis
y_range =   0:0.05:20;       % y-axis
[X, Y] = meshgrid(x_range, y_range);

% ---------- Initialize complex vector sum ----------
Re = zeros(size(X));
Im = zeros(size(X));

% ---------- Superpose N antennas (including phase progression) ----------
% Convention: the leftmost one is antenna 1, extra phase = (k-1)*Δ
for k = 1:N
    r = sqrt( (X - ant_x(k)).^2 + (Y - ant_y(k)).^2 );
    r(r==0) = eps;                              % Avoid division by zero
    Pden_k = Pt ./ (4*pi*r.^2);                 % Power density of a single antenna
    V_k = sqrt(Pden_k);                         % V = sqrt(Pden)
    p_geom = 2*pi*r/lambda;                     % Geometric phase
    p_extra = (k-1) * delta_rad;                % Phase difference between adjacent array elements
    p_k = p_geom + p_extra;                     % Total phase

    Re = Re + V_k .* cos(p_k);
    Im = Im + V_k .* sin(p_k);
end

% ---------- Combine and convert to dB ----------
Vabs = sqrt(Re.^2 + Im.^2);
Pden_total = Vabs.^2;            % Linear
PdB = 10*log10(Pden_total);      % dB(W/m^2)

% ---------- 2D colour map ----------
figure;
imagesc(x_range, y_range, PdB);
axis xy; axis equal; box on;
colormap(jet);
colorbar;
clim([-100 30]);  % Display range for easier comparison
title(sprintf('Phased Array Interference Pattern (Phase difference 60 degree) ', N, delta_deg));
xlabel('x (m)'); ylabel('y (m)');

% Mark antenna positions
hold on;
plot(ant_x, ant_y, 'wo', 'MarkerFaceColor','w', 'MarkerSize',6);