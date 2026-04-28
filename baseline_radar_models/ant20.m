clc; clear; close all;

% ---------- Parameters ----------
lambda = 0.012;          % Wavelength (m)
d = lambda/2;            % Spacing between adjacent antennas
N = 20;                  % Number of antennas
Pt = 1;                  % Transmit power of each antenna (W)

% Antenna coordinates: symmetric about the y-axis, with no element at the origin when the number is even
idx = (-(N-1)/2 : (N-1)/2);  % -9.5:9.5 (step size 1) -> total of N elements
ant_x = idx * d;             % x position
ant_y = zeros(size(ant_x));  % All y=0

% Calculate grid
x_range = -10:0.05:10;       % x-axis
y_range =   0:0.05:20;       % y-axis
[X, Y] = meshgrid(x_range, y_range);

% ---------- Initialize complex vector sum ----------
Re = zeros(size(X));
Im = zeros(size(X));

% ---------- Superpose N antennas ----------
for k = 1:N
    r = sqrt( (X - ant_x(k)).^2 + (Y - ant_y(k)).^2 );
    r(r==0) = eps;                              % Avoid division by zero
    Pden_k = Pt ./ (4*pi*r.^2);                 % Power density of a single antenna
    V_k = sqrt(Pden_k);                         % V = sqrt(Pden)
    p_k = 2*pi*r/lambda;                        % Phase 2πr/λ
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
clim([-100 30]);
title(sprintf('20-Antenna Interference Pattern (dB), N=%d', N));
xlabel('x (m)'); ylabel('y (m)');