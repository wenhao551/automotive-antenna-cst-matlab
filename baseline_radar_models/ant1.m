clc; clear; close all;

%-----------------------------
% Parameter settings
%-----------------------------
Pt = 1;                 % Transmit power (W)
x_range = -10:0.15:10;   % x-axis range
y_range = 0:0.15:20;     % y-axis range

%-----------------------------
% Generate grid coordinates
%-----------------------------
[X, Y] = meshgrid(x_range, y_range);

%-----------------------------
% Calculate distance r
% The antenna is located at the origin (0,0)
%-----------------------------
r = sqrt(X.^2 + Y.^2);
r(r == 0) = eps;        % Avoid division by zero

%-----------------------------
% Calculate power density (linear value)
% Pden = Pt / (4*pi*r^2)
%-----------------------------
Pden = Pt ./ (4 * pi * r.^2);

%-----------------------------
% Convert to dB scale
% PdB = 10*log10(Pden)
%-----------------------------
PdB = 10 * log10(Pden);

%-----------------------------
% Plot the 3D distribution (dB)
%-----------------------------
figure;
surf(X, Y, PdB);
shading interp;
colorbar;
title('Power Density Distribution (dB scale)');
xlabel('X axis (m)');
ylabel('Y axis (m)');
zlabel('Power Density P_{den} (dB, W/m^2)');
view(45, 30);

%-----------------------------
% Plot the contour map (dB)
%-----------------------------
figure;
contourf(X, Y, PdB, 30, 'LineColor', 'none');
colorbar;
title('Power Density Contour Map (dB scale)');
xlabel('X axis (m)');
ylabel('Y axis (m)');