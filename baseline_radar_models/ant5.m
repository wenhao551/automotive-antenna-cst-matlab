clc; clear; close all;

% ---------- Parameters ----------
Pt = 1;                  % Transmit power of each antenna (W)
lambda = 0.012;          % Wavelength (m)
d = lambda/2;            % Spacing between adjacent antennas
ant_x = [-2*d, -1*d, 0, +1*d, +2*d];  % x-positions of the five antennas (symmetric about the y-axis)
ant_y = zeros(size(ant_x));           % All antennas are at y=0

x_range = -10:0.15:10;   % x-axis
y_range =   0:0.15:20;   % y-axis
[X, Y] = meshgrid(x_range, y_range);

% ---------- Initialize complex vector sum ----------
Re = zeros(size(X));
Im = zeros(size(X));

% ---------- Superpose the fields from 5 antennas ----------
for k = 1:numel(ant_x)
    r = sqrt( (X - ant_x(k)).^2 + (Y - ant_y(k)).^2 );
    r(r==0) = eps;                               % Avoid division by zero
    Pden_k = Pt ./ (4*pi*r.^2);                  % Power density of a single antenna
    V_k = sqrt(Pden_k);                          % According to your formula: V = sqrt(Pden)
    p_k = 2*pi*r/lambda;                         % Phase
    Re = Re + V_k .* cos(p_k);
    Im = Im + V_k .* sin(p_k);
end

Vabs = sqrt(Re.^2 + Im.^2);          % Combined amplitude
Pden_total = Vabs.^2;                % Total power density (linear)
PdB = 10*log10(Pden_total);          % dB scale

% ---------- 2D colour map ----------
figure;
imagesc(x_range, y_range, PdB);
axis xy; axis equal; box on;
colormap(jet);
colorbar;

clim([-100 30]);

title('5-Antenna Interference Pattern (dB), N = 5');
xlabel('x (m)'); ylabel('y (m)');

% Mark antenna positions
hold on;
plot(ant_x, ant_y, 'wo', 'MarkerFaceColor','w', 'MarkerSize',8);
for k = 1:numel(ant_x)
    text(ant_x(k)+0.02, ant_y(k)-0.2, sprintf('A%d',k), 'Color','w');
end