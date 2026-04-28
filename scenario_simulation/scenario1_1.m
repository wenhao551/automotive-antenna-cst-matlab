clc; clear; close all;

%% =========================
% Radar / array parameters
% =========================
c = 3e8;                 % Speed of light (m/s)
f0 = 25e9;               % Radar centre frequency, 25 GHz

lambda0 = c / f0;        % Wavelength at the centre frequency, lambda = c / f
d = lambda0/2;           % Element spacing, set to lambda/2 to reduce obvious grating lobes
N = 20;                  % Number of array elements

Pt_total = 100;          % Total transmit power of the array (W)
Pt_elem  = Pt_total / N; % Transmit power assigned to each element

% n is used to generate the weighting window
n = 0:N-1;  
w = 0.54 - 0.46*cos(2*pi*n/(N-1));   % Hamming weighting

% The Hamming taper is used to reduce sidelobes.
% The trade-off is a slightly wider main lobe and a small peak loss.
w = w / sqrt(mean(w.^2));            % RMS normalization to keep the power scale similar

% mC is the centred element index relative to the array centre.
% For 20 elements, it is -9.5, -8.5, ..., 9.5.
mC = (0:N-1) - (N-1)/2;

%% =========================
% Ego vehicle
% =========================
L0 = 4.762;              % Ego vehicle length (m)
W0 = 1.847;              % Ego vehicle width (m)
xc0 = 0;                 
yc0 = 0;                 % Ego vehicle centre position

% Front and rear radar phase centres
xrF = xc0;  
yrF = yc0 + L0/2;

xrR = xc0;  
yrR = yc0 - L0/2;

% Front and rear array element coordinates
% Both arrays are arranged along the x direction.
ant_xF = xrF + mC*d;
ant_yF = yrF + zeros(size(mC));

ant_xR = xrR + mC*d;
ant_yR = yrR + zeros(size(mC));

% Ego vehicle rectangular outline for plotting
Erect.x = [-W0/2, W0/2, W0/2, -W0/2, -W0/2];
Erect.y = [-L0/2, -L0/2, L0/2, L0/2, -L0/2];

%% =========================
% Element pattern
% =========================
% Read the element radiation pattern exported from CST.
% Column 1: theta
% Column 2: phi
% Column 3: gain in dBi
T = readtable('CIR_Gain_Phi90.xlsx');

theta_all = T{:,1};
phi_all   = T{:,2};
gain_dBi  = T{:,3};

% Extract the phi = 90 degree cut with theta from 0 to 180 degrees.
mask = abs(phi_all - 90) < 1e-6 & theta_all >= 0 & theta_all <= 180;

theta_deg = theta_all(mask);

% Convert gain from dBi to linear scale.
G_lin = 10.^(gain_dBi(mask)/10);

% Sort by angle for stable interpolation.
[theta_deg, idx_sort] = sort(theta_deg);
G_lin = G_lin(idx_sort);

% Remove duplicated angle values to avoid interp1 errors.
[theta_deg, idx_unique] = unique(theta_deg, 'stable');
G_lin = G_lin(idx_unique);

% Element gain interpolation function.
% Input: local angle relative to the element boresight, in radians.
% Output: linear gain in that direction.
element_gain = @(a_rad) interp1(theta_deg, G_lin, abs(rad2deg(a_rad)), 'linear', 0);

%% =========================
% Multi-frequency averaging
% =========================
Nf = 11;                 % Number of frequency samples
B  = 400e6;              % Bandwidth, 400 MHz

% Generate frequency samples around the centre frequency.
f_vec      = f0 + linspace(-B/2, B/2, Nf);
lambda_vec = c ./ f_vec;
k_vec      = 2*pi ./ lambda_vec;   % Wavenumber, k = 2*pi / lambda

%% =========================
% Target parameters
% =========================
Lcar = L0;               % Target vehicle length
Wcar = W0;               % Target vehicle width
Lbik = 1.74;             % Simplified bicycle line length
step = 0.2;              % Scattering point sampling interval

%% =========================
% FRONT targets
% =========================
% Front vehicle 1
xcF1 = -7.5;
ycF1 = 12.0;

% mode = 'RB' means Right + Bottom edges are used as scattering points.
[XF1, YF1, RectF1] = make_rect_edges(xcF1, ycF1, Lcar, Wcar, step, 'RB');

% Front vehicle 2
xcF2 = 7.5;
ycF2 = 12.5;
[XF2, YF2, RectF2] = make_rect_edges(xcF2, ycF2, Lcar, Wcar, step, 'LB');

% Combine all front scattering points.
XptsF = [XF1 XF2];
YptsF = [YF1 YF2];

% typeF is used to distinguish different front targets.
% 1 = front vehicle 1
% 2 = front vehicle 2
typeF = [ones(1,numel(XF1)), 2*ones(1,numel(XF2))];

% Set the RCS of each scattering point to 1.
sigmaF = ones(size(XptsF));

% No extra shadowing is applied to the front targets here.
shadowF = false(size(XptsF));

%% =========================
% REAR targets
% =========================
% Rear vehicle
[XR1, YR1, RectR1, InfoR1] = make_rect_edges(-6.5, -12.0, Lcar, Wcar, step, 'RT');

% Rear-left bicycle
[XR2, YR2, BikeL2] = make_bicycle_line(-8.6, -15.0, Lbik, step, 'vertical');

% Rear-right bicycle
[XR3, YR3, BikeL3] = make_bicycle_line( 6.5, -14.6, Lbik, step, 'vertical');

% Combine all rear scattering points.
XptsR = [XR1 XR2 XR3];
YptsR = [YR1 YR2 YR3];

% 1 = rear vehicle
% 2 = rear-left bicycle
% 3 = rear-right bicycle
typeR = [ones(1,numel(XR1)), 2*ones(1,numel(XR2)), 3*ones(1,numel(XR3))];

sigmaR = ones(size(XptsR));

% Initialize rear shadow mask.
shadowR = false(size(XptsR));

% Only check whether the two bicycles are blocked by the rear vehicle.
shadowR(typeR ~= 1) = mark_shadow_by_blocker( ...
    xrR, yrR, XptsR(typeR ~= 1), YptsR(typeR ~= 1), InfoR1);

% Print the number of rear shadowed points for checking.
fprintf('Rear shadowed points: %d / %d\n', sum(shadowR), numel(XptsR));

%% =========================
% Random phase
% =========================
% Assign a random initial phase to each scattering point.
rng(1);                         % Fix random seed for repeatable results
phi_rand_F = 2*pi*rand(size(XptsF));
phi_rand_R = 2*pi*rand(size(XptsR));

%% =========================
% Scan angles
% =========================
Nd = 1081;                     % Number of scan angle samples

% Front array scans from 0 to 180 degrees.
theta_front = linspace(0, pi, Nd);

% Rear array scans from 180 to 360 degrees.
theta_rear  = linspace(pi, 2*pi, Nd);

%% =========================
% Echo power
% =========================
% Calculate echo power for each front scattering point at each scan angle.
Pr_each_F = compute_Pr_each( ...
    XptsF, YptsF, sigmaF, shadowF, xrF, yrF, ant_xF, ant_yF, ...
    theta_front, d, N, Pt_total, Pt_elem, w, phi_rand_F, ...
    k_vec, lambda_vec, element_gain, +1);

% Calculate echo power for each rear scattering point at each scan angle.
% facing_sign = -1 means the rear array faces the opposite direction.
Pr_each_R = compute_Pr_each( ...
    XptsR, YptsR, sigmaR, shadowR, xrR, yrR, ant_xR, ant_yR, ...
    theta_rear, d, N, Pt_total, Pt_elem, w, phi_rand_R, ...
    k_vec, lambda_vec, element_gain, -1);

%% =========================
% Polar accumulation
% =========================
r_step = 0.2;
% Range resolution of the polar grid

x_min = -15; x_max = 15;
y_min = -20; y_max = 20;
% Final display range in the x-y plane

corners = [x_min y_min;
           x_min y_max;
           x_max y_min;
           x_max y_max];
% Corner coordinates of the display area

[XplotF, YplotF, PplotF, ~] = build_polar_map( ...
    XptsF, YptsF, Pr_each_F, xrF, yrF, theta_front, ...
    r_step, x_min, x_max, y_min, y_max, corners);

[XplotR, YplotR, PplotR, ~] = build_polar_map( ...
    XptsR, YptsR, Pr_each_R, xrR, yrR, theta_rear, ...
    r_step, x_min, x_max, y_min, y_max, corners);

global_max = max([PplotF(:); PplotR(:)]);
% Maximum displayed power from both front and rear arrays

dyn_range = 40;
% Display dynamic range, in dB

%% =========================
% Plot
% =========================
figure;
clf;
hold on;
axis equal;      % Keep the same scale for x and y
grid on;         % Enable grid
box on;          % Show plot box

% Set figure background to white.
set(gcf, 'Color', 'w');

% Set colormap.
cmap = jet(256);
colormap(cmap);

% Get current axes.
ax = gca;

% Set the axes background to the lowest colour in the colormap.
ax.Color = cmap(1,:);

% Draw x-axis reference line.
plot([x_min x_max], [0 0], 'k-', 'LineWidth', 0.8);

% Draw y-axis reference line.
plot([0 0], [y_min y_max], 'k-', 'LineWidth', 0.8);

% Set display limits.
xlim([x_min x_max]);
ylim([y_min y_max]);

% Front power distribution.
scatter(XplotF, YplotF, 28, PplotF, 's', 'filled', 'MarkerEdgeColor', 'none');

% Rear power distribution.
scatter(XplotR, YplotR, 28, PplotR, 's', 'filled', 'MarkerEdgeColor', 'none');

% Overlay target outlines, scattering points, array elements, and ego vehicle outline.
draw_scene(RectF1, RectF2, RectR1, BikeL2, BikeL3, ...
           XF1, YF1, XF2, YF2, XR1, YR1, XR2, YR2, XR3, YR3, ...
           XptsR(shadowR), YptsR(shadowR), ...
           ant_xF, ant_yF, ant_xR, ant_yR, ...
           Erect, xrF, yrF, xrR, yrR);

caxis([global_max-dyn_range, global_max]);
% Set colour scale.

cb = colorbar;
ylabel(cb, 'Received power, P_r (dBW)');

title('Scenario 1');
xlabel('Lateral position, x (m)');
ylabel('Longitudinal position, y (m)');

hold off;

%% =========================
% Functions
% =========================

function [X, Y, rect, info] = make_rect_edges(xc, yc, L, W, step, mode)
    % Generate discrete scattering points on two selected edges of a rectangle.
    %
    % 'LB' = Left + Bottom
    % 'RB' = Right + Bottom
    % 'LT' = Left + Top
    % 'RT' = Right + Top

    build = @(a,b) linspace(a, b, round((b-a)/step)+1);

    xL = xc - W/2;
    xR = xc + W/2;
    yB = yc - L/2;
    yT = yc + L/2;

    switch upper(mode)
        case 'LB'
            xx = [xL*ones(size(build(yB,yT))), build(xL,xR)];
            yy = [build(yB,yT), yB*ones(size(build(xL,xR)))];

        case 'RB'
            xx = [xR*ones(size(build(yB,yT))), build(xL,xR)];
            yy = [build(yB,yT), yB*ones(size(build(xL,xR)))];

        case 'LT'
            xx = [xL*ones(size(build(yB,yT))), build(xL,xR)];
            yy = [build(yB,yT), yT*ones(size(build(xL,xR)))];

        case 'RT'
            xx = [xR*ones(size(build(yB,yT))), build(xL,xR)];
            yy = [build(yB,yT), yT*ones(size(build(xL,xR)))];

        otherwise
            error('Unknown mode');
    end

    XY = unique([xx(:), yy(:)], 'rows', 'stable');

    X = XY(:,1).';
    Y = XY(:,2).';

    rect.x = [xL xR xR xL xL];
    rect.y = [yB yB yT yT yB];

    info.x_left  = xL;
    info.x_right = xR;
    info.y_bot   = yB;
    info.y_top   = yT;
end

function [X, Y, bikeLine] = make_bicycle_line(xc, yc, Lbik, step, orientation)
    % Generate a simplified bicycle line target.

    s = -Lbik/2 : step : Lbik/2;

    if abs(s(end)-Lbik/2) > 1e-9
        s(end+1) = Lbik/2;
    end

    switch lower(orientation)
        case 'vertical'
            X = xc + zeros(size(s));
            Y = yc + s;

        case 'horizontal'
            X = xc + s;
            Y = yc + zeros(size(s));

        otherwise
            error('Unknown bicycle orientation');
    end

    bikeLine.x = X;
    bikeLine.y = Y;
end

function shadow = mark_shadow_by_blocker(xr, yr, X, Y, info)
    % Check whether target points are blocked by the blocker vehicle.

    tol = 1e-9;

    p1 = [info.x_right, info.y_bot];
    p2 = [info.x_left,  info.y_top];
    v = p2 - p1;

    radar_vec  = [xr - p1(1), yr - p1(2)];
    radar_side = v(1)*radar_vec(2) - v(2)*radar_vec(1);

    th1 = atan2(p1(2)-yr, p1(1)-xr);
    th2 = atan2(p2(2)-yr, p2(1)-xr);

    if abs(th1-th2) > pi
        if th1 < th2
            th1 = th1 + 2*pi;
        else
            th2 = th2 + 2*pi;
        end
    end

    th_min = min(th1, th2);
    th_max = max(th1, th2);

    shadow = false(size(X));

    for k = 1:numel(X)
        th = atan2(Y(k)-yr, X(k)-xr);

        if th < th_min-pi
            th = th + 2*pi;
        elseif th > th_max+pi
            th = th - 2*pi;
        end

        if th < th_min || th > th_max
            continue;
        end

        w    = [X(k)-p1(1), Y(k)-p1(2)];
        side = v(1)*w(2) - v(2)*w(1);

        if (side * radar_side < 0) || (abs(side) < tol)
            shadow(k) = true;
        end
    end
end

function Pr_each = compute_Pr_each( ...
    Xpts, Ypts, sigma, is_shadow, xr, yr, ant_x, ant_y, ...
    theta_scan, d, N, Pt_total, Pt_elem, w, phi_rand, ...
    k_vec, lambda_vec, element_gain, facing_sign)

    Nd = numel(theta_scan);
    Nf = numel(k_vec);
    Nt = numel(Xpts);

    mC = (0:N-1) - (N-1)/2;
    Rpts = hypot(Xpts - xr, Ypts - yr);
    Piso = Pt_elem ./ (4*pi*Rpts.^2);

    Pr_each = zeros(Nt, Nd);

    for fi = 1:Nf
        kf     = k_vec(fi);
        lambda = lambda_vec(fi);
        dphi = kf * d * cos(theta_scan);

        for i = 1:Nd
            phi_feed = mC * dphi(i);

            for t = 1:Nt
                if is_shadow(t)
                    continue;
                end

                dx = Xpts(t) - ant_x;
                dy = Ypts(t) - ant_y;
                rk = hypot(dx, dy);

                a_k = atan2(dx, facing_sign * dy);
                G_elem = element_gain(a_k);
                Pden_k = G_elem .* Pt_elem ./ (4*pi*rk.^2);
                Vk = w .* sqrt(Pden_k);

                phi = kf*rk + phi_feed + phi_rand(t);
                E_sum = sum(Vk .* exp(1j*phi));
                P_array = abs(E_sum)^2;
                G = P_array / (N * Piso(t));

                Pr_each(t,i) = Pr_each(t,i) + ...
                    Pt_total * G^2 * lambda^2 * sigma(t) / ((4*pi)^3 * Rpts(t)^4);
            end
        end
    end

    Pr_each = Pr_each / Nf;
end

function [Xplot, Yplot, Prplot, maxVal] = build_polar_map( ...
    Xpts, Ypts, Pr_each, xr, yr, theta_scan, ...
    r_step_cell, x_min, x_max, y_min, y_max, corners)

    Nd = numel(theta_scan);
    Rmax = max(hypot(corners(:,1)-xr, corners(:,2)-yr));
    Nr = floor(Rmax / r_step_cell);
    r_centers = ((1:Nr) - 0.5) * r_step_cell;

    Rpts = hypot(Xpts - xr, Ypts - yr);
    r_idx = floor(Rpts / r_step_cell) + 1;
    r_idx = max(1, min(Nr, r_idx));

    theta_actual = atan2(Ypts - yr, Xpts - xr);
    theta_scan_min = min(theta_scan);
    theta_scan_max = max(theta_scan);

    Pr_polar = zeros(Nr, Nd);

    for t = 1:numel(Xpts)
        theta_t = theta_actual(t);

        if theta_t < theta_scan_min
            theta_t = theta_t + 2*pi;
        end
        if theta_t > theta_scan_max && theta_t < theta_scan_min + 2*pi
            theta_t = theta_t - 2*pi;
        end

        if theta_t < theta_scan_min || theta_t > theta_scan_max
            continue;
        end

        [~, idx_angle] = min(abs(theta_scan - theta_t));

        angle_window_deg = 10;
        angle_window_idx = round(angle_window_deg / rad2deg(theta_scan(2)-theta_scan(1)));

        for i = max(1, idx_angle-angle_window_idx) : min(Nd, idx_angle+angle_window_idx)
            Pr_polar(r_idx(t), i) = Pr_polar(r_idx(t), i) + Pr_each(t, i);
        end
    end

    Pr_polar_dBW = 10*log10(Pr_polar + realmin);

    [ThetaGrid, RGrid] = meshgrid(theta_scan, r_centers);
    Xgrid = xr + RGrid .* cos(ThetaGrid);
    Ygrid = yr + RGrid .* sin(ThetaGrid);

    mask = (Xgrid >= x_min) & (Xgrid <= x_max) & ...
           (Ygrid >= y_min) & (Ygrid <= y_max);

    Xplot  = Xgrid(mask);
    Yplot  = Ygrid(mask);
    Prplot = Pr_polar_dBW(mask);

    maxVal = max(Prplot(:));
end

function draw_scene(RectF1, RectF2, RectR1, BikeL2, BikeL3, ...
                    XF1, YF1, XF2, YF2, XR1, YR1, XR2, YR2, XR3, YR3, ...
                    Xsh, Ysh, ant_xF, ant_yF, ant_xR, ant_yR, ...
                    Erect, xrF, yrF, xrR, yrR)

    plot(RectF1.x, RectF1.y, 'k-', 'LineWidth', 2.0);
    plot(RectF2.x, RectF2.y, 'k-', 'LineWidth', 2.0);
    plot(RectR1.x, RectR1.y, 'k-', 'LineWidth', 2.0);

    plot(BikeL2.x, BikeL2.y, 'b-', 'LineWidth', 2.4);
    plot(BikeL3.x, BikeL3.y, 'b-', 'LineWidth', 2.4);

    plot(XF1, YF1, 'ko', 'MarkerFaceColor','y', 'MarkerSize',5);
    plot(XF2, YF2, 'ks', 'MarkerFaceColor','y', 'MarkerSize',5);

    plot(XR1, YR1, 'ko', 'MarkerFaceColor','c', 'MarkerSize',5);
    plot(XR2, YR2, 'bd', 'MarkerFaceColor','c', 'MarkerSize',5);
    plot(XR3, YR3, 'bd', 'MarkerFaceColor','c', 'MarkerSize',5);

    plot(Xsh, Ysh, 'rx', 'MarkerSize', 8, 'LineWidth', 1.2);

    plot(ant_xF, ant_yF, 'm.', 'MarkerSize', 6);
    plot(ant_xR, ant_yR, 'm.', 'MarkerSize', 6);

    plot(Erect.x, Erect.y, 'm-', 'LineWidth', 2.5);
    plot(xrF, yrF, 'mp', 'MarkerFaceColor','m', 'MarkerSize',12);
    plot(xrR, yrR, 'mp', 'MarkerFaceColor','m', 'MarkerSize',12);
end