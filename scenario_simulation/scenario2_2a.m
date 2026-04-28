clc; clear; close all;
% clc        : Clear previous text in the command window
% clear      : Clear existing variables in the workspace
% close all  : Close all open figure windows

%% =========================
% Radar / array parameters
% =========================
c = 3e8;
% Speed of light in m/s

f0 = 25e9;
% Centre operating frequency of the radar / antenna system, in Hz

lambda0 = c / f0;
% Wavelength at the centre frequency

d = lambda0/2;
% Element spacing, set to half wavelength

N = 20;
% Number of antenna elements

Pt_total = 100;
% Total transmit power of the array

Pt_elem  = Pt_total / N;
% Transmit power assigned to each element

n = 0:N-1;
% Element index

w = 0.54 - 0.46*cos(2*pi*n/(N-1));
% Hamming taper used to reduce sidelobes

w = w / sqrt(mean(w.^2));
% Normalize the taper weights

mC = (0:N-1) - (N-1)/2;
% Centred element index relative to the array centre

%% =========================
% Ego vehicle
% =========================
L0 = 4.762;
% Ego vehicle length, in metres

W0 = 1.847;
% Ego vehicle width, in metres

xc0 = 0; yc0 = 0;
% Ego vehicle centre position

% Front and rear radar phase centres
xrF = xc0;  yrF = yc0 + L0/2;
% Front radar phase centre

xrR = xc0;  yrR = yc0 - L0/2;
% Rear radar phase centre

% Front and rear array element coordinates
ant_xF = xrF + mC*d; ant_yF = yrF + zeros(size(mC));
% Front array arranged along the x direction

ant_xR = xrR + mC*d; ant_yR = yrR + zeros(size(mC));
% Rear array arranged along the x direction

% Ego vehicle outline
xL0 = xc0 - W0/2;
% Left boundary of the ego vehicle

xR0 = xc0 + W0/2;
% Right boundary of the ego vehicle

yB0 = yc0 - L0/2;
% Rear boundary of the ego vehicle

yT0 = yc0 + L0/2;
% Front boundary of the ego vehicle

Erect.x = [xL0 xR0 xR0 xL0 xL0];
Erect.y = [yB0 yB0 yT0 yT0 yB0];
% Rectangular outline of the ego vehicle

%% =========================
% Element pattern
% =========================
T = readtable('CIR_Gain_Phi90.xlsx');
% Read the CST-exported antenna gain data

theta_all = T{:,1};
% Column 1: theta angle

phi_all   = T{:,2};
% Column 2: phi angle

gain_dBi  = T{:,3};
% Column 3: gain in dBi

mask = abs(phi_all - 90) < 1e-6 & theta_all >= 0 & theta_all <= 180;
% Extract the phi = 90 degree cut with theta from 0 to 180 degrees

theta_deg = theta_all(mask);
% Selected theta values

G_lin = 10.^(gain_dBi(mask)/10);
% Convert gain from dBi to linear scale

[theta_deg, idx] = sort(theta_deg);
% Sort angles for interpolation

G_lin = G_lin(idx);
% Reorder gain values using the same order

[theta_deg, idx2] = unique(theta_deg, 'stable');
% Remove duplicated angle values

G_lin = G_lin(idx2);
% Keep the corresponding unique gain values

element_gain = @(a_rad) interp1(theta_deg, G_lin, abs(rad2deg(a_rad)), 'linear', 0);
% Interpolation function for element gain

%% =========================
% Multi-frequency averaging
% =========================
Nf = 11;
% Number of frequency samples

B  = 400e6;
% System bandwidth, 400 MHz

f_vec      = f0 + linspace(-B/2, B/2, Nf);
% Frequency samples around the centre frequency

lambda_vec = c ./ f_vec;
% Wavelength at each frequency sample

k_vec      = 2*pi ./ lambda_vec;
% Wavenumber at each frequency sample

%% =========================
% Target parameters
% =========================
Lcar = L0;
% Target vehicle length

Wcar = W0;
% Target vehicle width

Lbik = 1.74;
% Simplified bicycle length

step = 0.2;
% Sampling interval along target edges

%% =========================
% FRONT targets
% =========================
% Front target: one upper-left vehicle
[XF1, YF1, RectF1] = make_rect_edges(-7.5, 12.0, Lcar, Wcar, step, 'RB');
% Generate scattering points on the right and bottom edges

XptsF = XF1;
YptsF = YF1;
% Front target scattering points

typeF = ones(1, numel(XF1));
% Target label for the front vehicle

sigmaF = ones(size(XptsF));
% RCS of each scattering point

shadowF = false(size(XptsF));
% Initial front shadow mask

%% =========================
% REAR / SIDE targets
% =========================
% Rear target: one lower-left vehicle
[XR1, YR1, RectR1, InfoR1] = make_rect_edges(-6.5, -12.0, Lcar, Wcar, step, 'RT');
% Generate scattering points on the right and top edges

% Right-side target: one bicycle
[XR2, YR2, BikeR] = make_bicycle_line(6.5, 0.0, Lbik, step, 'vertical');
% Generate a vertical line target for the bicycle

XptsR = [XR1 XR2];
YptsR = [YR1 YR2];
% Combined rear / side scattering points

% 1 = rear vehicle, 2 = right-side bicycle
typeR = [ones(1,numel(XR1)), 2*ones(1,numel(XR2))];
% Target labels

sigmaR = ones(size(XptsR));
% RCS of each scattering point

shadowR = false(size(XptsR));
% Initial rear / side shadow mask

%% =========================
% External shadow for rear scene
% =========================
shadowR(typeR ~= 1) = mark_shadow_by_blocker( ...
    xrR, yrR, XptsR(typeR ~= 1), YptsR(typeR ~= 1), InfoR1);
% Check whether the bicycle is blocked by the rear vehicle

%% =========================
% Ego self-shadow
% =========================
selfShadowF = false(size(XptsF));
% Front target self-shadow mask caused by the ego vehicle

for k = 1:numel(XptsF)
    % Check each front target point
    selfShadowF(k) = los_blocked_by_ego( ...
        xrF, yrF, XptsF(k), YptsF(k), xL0, xR0, yB0, yT0);
end

selfShadowR = false(size(XptsR));
% Rear / side target self-shadow mask caused by the ego vehicle

for k = 1:numel(XptsR)
    selfShadowR(k) = los_blocked_by_ego( ...
        xrR, yrR, XptsR(k), YptsR(k), xL0, xR0, yB0, yT0);
end

% Combine shadow masks for Pr calculation
shadowF = shadowF | selfShadowF;
% Final front shadow mask

shadowR = shadowR | selfShadowR;
% Final rear / side shadow mask

fprintf('Front self-shadowed points: %d / %d\n', sum(selfShadowF), numel(XptsF));
% Print number of front self-shadowed points

fprintf('Rear shadowed points: %d / %d\n', sum(shadowR), numel(XptsR));
% Print number of rear / side shadowed points

fprintf('Right-side bicycle shadowed points (used in Pr): %d / %d\n', ...
    sum(shadowR(typeR==2)), sum(typeR==2));
% Print number of shadowed bicycle points used in Pr calculation

%% =========================
% Plot-only shadow selection
% =========================
% Bicycle shadowing is used in Pr calculation, but not marked with red crosses in the plot
plotShadowR = shadowR & (typeR == 1);
% Only show rear-vehicle shadow points in the final plot

%% =========================
% Random phase
% =========================
rng(1);
% Fix random seed for repeatable results

phi_rand_F = 2*pi*rand(size(XptsF));
% Random initial phase for front scattering points

phi_rand_R = 2*pi*rand(size(XptsR));
% Random initial phase for rear / side scattering points

%% =========================
% Scan angles
% =========================
Nd = 1081;
% Number of scan angle samples

theta_front = linspace(0, pi, Nd);
% Front scan range: 0 to 180 degrees

theta_rear  = linspace(pi, 2*pi, Nd);
% Rear scan range: 180 to 360 degrees

%% =========================
% Echo power
% =========================
Pr_each_F = compute_Pr_each( ...
    XptsF, YptsF, sigmaF, shadowF, xrF, yrF, ant_xF, ant_yF, ...
    theta_front, d, N, Pt_total, Pt_elem, w, phi_rand_F, ...
    k_vec, lambda_vec, element_gain, +1);
% Calculate received power for the front target points

Pr_each_R = compute_Pr_each( ...
    XptsR, YptsR, sigmaR, shadowR, xrR, yrR, ant_xR, ant_yR, ...
    theta_rear, d, N, Pt_total, Pt_elem, w, phi_rand_R, ...
    k_vec, lambda_vec, element_gain, -1);
% Calculate received power for the rear / side target points

%% =========================
% Polar accumulation
% =========================
r_step = 0.2;
% Range resolution of the polar grid

x_min = -15; x_max = 15;
y_min = -20; y_max = 20;
% Display range in the x-y plane

corners = [x_min y_min; x_min y_max; x_max y_min; x_max y_max];
% Corner coordinates of the display area

[XplotF, YplotF, PplotF, ~] = build_polar_map( ...
    XptsF, YptsF, Pr_each_F, xrF, yrF, theta_front, ...
    r_step, x_min, x_max, y_min, y_max, corners);
% Convert front polar data to x-y plot points

[XplotR, YplotR, PplotR, ~] = build_polar_map( ...
    XptsR, YptsR, Pr_each_R, xrR, yrR, theta_rear, ...
    r_step, x_min, x_max, y_min, y_max, corners);
% Convert rear polar data to x-y plot points

global_max = max([PplotF(:); PplotR(:)]);
% Maximum displayed power value

dyn_range = 40;
% Display dynamic range in dB

%% =========================
% Plot
% =========================
figure;
% Create a new figure

hold on;
% Allow multiple objects on the same plot

axis equal;
% Keep equal scale on x and y axes

grid on;
% Show grid

box on;
% Show plot box

set(gcf,'Color','w');
% Set figure background to white

cmap = jet(256);
% Generate jet colormap

colormap(cmap);
% Apply colormap

ax = gca;
% Get current axes

ax.Color = cmap(1,:);
% Set axes background to the lowest colormap colour

plot([x_min x_max],[0 0],'k-','LineWidth',0.8);
% Draw x-axis reference line

plot([0 0],[y_min y_max],'k-','LineWidth',0.8);
% Draw y-axis reference line

xlim([x_min x_max]);
% Set x-axis range

ylim([y_min y_max]);
% Set y-axis range

scatter(XplotF, YplotF, 28, PplotF, 's', 'filled', 'MarkerEdgeColor','none');
% Plot front power distribution

scatter(XplotR, YplotR, 28, PplotR, 's', 'filled', 'MarkerEdgeColor','none');
% Plot rear power distribution

draw_scene(RectF1, RectR1, BikeR, ...
           XF1, YF1, XR1, YR1, XR2, YR2, ...
           XptsF(selfShadowF), YptsF(selfShadowF), ...
           XptsR(plotShadowR), YptsR(plotShadowR), ...
           ant_xF, ant_yF, ant_xR, ant_yR, ...
           Erect, xrF, yrF, xrR, yrR);
% Overlay target geometry, scattering points, shadow points, arrays, and ego vehicle

caxis([global_max-dyn_range, global_max]);
% Set colour scale

cb = colorbar;
% Add colourbar
title('Scenario 2a');
xlabel('Lateral position, x (m)');
ylabel('Longitudinal position, y (m)');
ylabel(cb, 'Received power, P_r (dBW)');


% Label the colourbar

% Plot title and axis labels

hold off;
% Stop adding objects to the current plot

%% =========================
% Functions
% =========================

function [X, Y, rect, info] = make_rect_edges(xc, yc, L, W, step, mode)
% Generate sampled edge points for a rectangular target.
%
% Inputs:
% xc, yc : rectangle centre
% L, W   : length and width
% step   : edge sampling interval
% mode   : selected scattering edges
%
% Outputs:
% X, Y   : sampled edge point coordinates
% rect   : complete rectangle outline for plotting
% info   : rectangle boundary information

    build = @(a,b) linspace(a, b, round((b-a)/step)+1);
    % Generate points between a and b

    xL = xc - W/2; xR = xc + W/2;
    % Left and right boundaries

    yB = yc - L/2; yT = yc + L/2;
    % Bottom and top boundaries

    switch upper(mode)
        % Select two edges according to mode

        case 'LB'
            xx = [xL*ones(size(build(yB,yT))), build(xL,xR)];
            yy = [build(yB,yT), yB*ones(size(build(xL,xR)))];
            % Left + bottom edges

        case 'RB'
            xx = [xR*ones(size(build(yB,yT))), build(xL,xR)];
            yy = [build(yB,yT), yB*ones(size(build(xL,xR)))];
            % Right + bottom edges

        case 'LT'
            xx = [xL*ones(size(build(yB,yT))), build(xL,xR)];
            yy = [build(yB,yT), yT*ones(size(build(xL,xR)))];
            % Left + top edges

        case 'RT'
            xx = [xR*ones(size(build(yB,yT))), build(xL,xR)];
            yy = [build(yB,yT), yT*ones(size(build(xL,xR)))];
            % Right + top edges

        otherwise
            error('Unknown mode');
            % Report unknown edge mode
    end

    XY = unique([xx(:), yy(:)], 'rows', 'stable');
    % Remove duplicated corner points

    X = XY(:,1).';
    Y = XY(:,2).';
    % Output x and y as row vectors

    rect.x = [xL xR xR xL xL];
    rect.y = [yB yB yT yT yB];
    % Complete rectangle outline

    info.x_left = xL;
    info.x_right = xR;
    info.y_bot = yB;
    info.y_top = yT;
    % Store boundary information
end

function [X, Y, bikeLine] = make_bicycle_line(xc, yc, Lbik, step, orientation)
% Generate a simplified bicycle line target.
%
% Inputs:
% xc, yc      : bicycle centre position
% Lbik        : bicycle length
% step        : sampling interval
% orientation : horizontal or vertical
%
% Outputs:
% X, Y        : sampled points on the bicycle line
% bikeLine    : line data for plotting

    s = -Lbik/2:step:Lbik/2;
    % Generate points along the bicycle length

    if abs(s(end)-Lbik/2) > 1e-9
        s(end+1) = Lbik/2;
        % Add the end point if needed
    end

    switch lower(orientation)
        case 'vertical'
            X = xc + zeros(size(s));
            Y = yc + s;
            % Vertical line

        case 'horizontal'
            X = xc + s;
            Y = yc + zeros(size(s));
            % Horizontal line

        otherwise
            error('Unknown bicycle orientation');
    end

    bikeLine.x = X;
    bikeLine.y = Y;
    % Store line data for plotting
end

function shadow = mark_shadow_by_blocker(xr, yr, X, Y, info)
% Check whether target points are blocked by a rectangular blocker.

    tol = 1e-9;
    % Numerical tolerance

    p1 = [info.x_right, info.y_bot];
    p2 = [info.x_left,  info.y_top];
    % Two diagonal points of the blocker

    v  = p2 - p1;
    % Diagonal direction vector

    radar_vec = [xr - p1(1), yr - p1(2)];
    % Vector from p1 to radar

    radar_side = v(1)*radar_vec(2) - v(2)*radar_vec(1);
    % Side of the diagonal where the radar is located

    th1 = atan2(p1(2)-yr, p1(1)-xr);
    th2 = atan2(p2(2)-yr, p2(1)-xr);
    % Angular limits of the blocker

    if abs(th1-th2) > pi
        % Handle angle wrapping
        if th1 < th2
            th1 = th1 + 2*pi;
        else
            th2 = th2 + 2*pi;
        end
    end

    th_min = min(th1, th2);
    th_max = max(th1, th2);
    % Shadow angular range

    shadow = false(size(X));
    % Initialise output mask

    for k = 1:numel(X)
        % Check each target point
        th = atan2(Y(k)-yr, X(k)-xr);
        % Target angle relative to radar

        if th < th_min-pi
            th = th + 2*pi;
        elseif th > th_max+pi
            th = th - 2*pi;
        end
        % Adjust angle range if needed

        if th < th_min || th > th_max
            continue;
            % Skip points outside the blocker angular sector
        end

        w = [X(k)-p1(1), Y(k)-p1(2)];
        % Vector from p1 to target point

        side = v(1)*w(2) - v(2)*w(1);
        % Side of the diagonal where the target point is located

        if (side * radar_side < 0) || (abs(side) < tol)
            shadow(k) = true;
            % Mark the point as shadowed
        end
    end
end

function blocked = los_blocked_by_ego(xr, yr, xt, yt, xL, xR, yB, yT)
% Check whether the line of sight from radar to target crosses the ego vehicle.

    tol = 1e-9;
    blocked = false;
    % Default: not blocked

    dx = xt - xr;
    dy = yt - yr;
    % Direction vector from radar to target

    if abs(dx) > tol
        % Check intersections with left and right boundaries

        t1 = (xL - xr)/dx;
        y1 = yr + t1*dy;
        % Intersection with left boundary

        if (t1 > tol) && (t1 < 1-tol) && (y1 > yB+tol) && (y1 < yT-tol)
            blocked = true;
            return;
            % Line of sight passes through the ego vehicle
        end

        t2 = (xR - xr)/dx;
        y2 = yr + t2*dy;
        % Intersection with right boundary

        if (t2 > tol) && (t2 < 1-tol) && (y2 > yB+tol) && (y2 < yT-tol)
            blocked = true;
            return;
        end
    end

    if abs(dy) > tol
        % Check intersections with bottom and top boundaries

        t3 = (yB - yr)/dy;
        x3 = xr + t3*dx;
        % Intersection with bottom boundary

        if (t3 > tol) && (t3 < 1-tol) && (x3 > xL+tol) && (x3 < xR-tol)
            blocked = true;
            return;
        end

        t4 = (yT - yr)/dy;
        x4 = xr + t4*dx;
        % Intersection with top boundary

        if (t4 > tol) && (t4 < 1-tol) && (x4 > xL+tol) && (x4 < xR-tol)
            blocked = true;
            return;
        end
    end
end

function Pr_each = compute_Pr_each( ...
    Xpts, Ypts, sigma, is_shadow, xr, yr, ant_x, ant_y, ...
    theta_scan, d, N, Pt_total, Pt_elem, w, phi_rand, ...
    k_vec, lambda_vec, element_gain, facing_sign)
% Calculate received power Pr for each target point and each scan angle.

    Nd = numel(theta_scan);
    % Number of scan angles

    Nf = numel(k_vec);
    % Number of frequency samples

    Nt = numel(Xpts);
    % Number of target points

    mC = (0:N-1) - (N-1)/2;
    % Centred element index

    Rpts = hypot(Xpts - xr, Ypts - yr);
    % Distance from radar phase centre to each target point

    Piso = Pt_elem ./ (4*pi*Rpts.^2);
    % Isotropic single-element reference power density

    Pr_each = zeros(Nt, Nd);
    % Output matrix

    for fi = 1:Nf
        % Loop over frequency samples

        kf = k_vec(fi);
        % Current wavenumber

        lambda = lambda_vec(fi);
        % Current wavelength

        dphi = kf * d * cos(theta_scan);
        % Progressive feed phase for beam scanning

        for i = 1:Nd
            % Loop over scan angles

            phi_feed = mC * dphi(i);
            % Feed phase of each element

            for t = 1:Nt
                % Loop over target points

                if is_shadow(t), continue; end
                % Skip shadowed points

                dx = Xpts(t) - ant_x;
                dy = Ypts(t) - ant_y;
                % Relative position from each element to the target point

                rk = hypot(dx, dy);
                % Distance from each element to the target point

                a_k = atan2(dx, facing_sign * dy);
                % Local angle relative to the array facing direction

                G_elem = element_gain(a_k);
                % Element gain in this direction

                Pden_k = G_elem .* Pt_elem ./ (4*pi*rk.^2);
                % Power density from each element at the target point

                Vk = w .* sqrt(Pden_k);
                % Convert power density to field amplitude and apply taper

                phi = kf*rk + phi_feed + phi_rand(t);
                % Total phase

                E_sum = sum(Vk .* exp(1j*phi));
                % Coherent field summation

                P_array = abs(E_sum)^2;
                % Array power quantity

                G = P_array / (N * Piso(t));
                % Equivalent array gain

                Pr_each(t,i) = Pr_each(t,i) + ...
                    Pt_total * G^2 * lambda^2 * sigma(t) / ((4*pi)^3 * Rpts(t)^4);
                % Monostatic radar equation
            end
        end
    end

    Pr_each = Pr_each / Nf;
    % Average over all frequency samples
end

function [Xplot, Yplot, Prplot, maxVal] = build_polar_map( ...
    Xpts, Ypts, Pr_each, xr, yr, theta_scan, ...
    r_step_cell, x_min, x_max, y_min, y_max, corners)
% Accumulate target-point power into a polar range-angle map,
% then convert it back to x-y coordinates for plotting.

    Nd = numel(theta_scan);
    % Number of scan angles

    Rmax = max(hypot(corners(:,1)-xr, corners(:,2)-yr));
    % Maximum display distance from the radar

    Nr = floor(Rmax / r_step_cell);
    % Number of range cells

    r_centers = ((1:Nr) - 0.5) * r_step_cell;
    % Range-cell centre positions

    Rpts = hypot(Xpts - xr, Ypts - yr);
    % Distance from each target point to radar

    r_idx = floor(Rpts / r_step_cell) + 1;
    % Range-cell index for each target point

    r_idx = max(1, min(Nr, r_idx));
    % Keep indices within valid range

    theta_actual = atan2(Ypts - yr, Xpts - xr);
    % Actual target angle relative to radar

    theta_scan_min = min(theta_scan);
    theta_scan_max = max(theta_scan);
    % Scan angle limits

    Pr_polar = zeros(Nr, Nd);
    % Polar power map

    for t = 1:numel(Xpts)
        theta_t = theta_actual(t);
        % Current target angle

        if theta_t < theta_scan_min
            theta_t = theta_t + 2*pi;
        end
        if theta_t > theta_scan_max && theta_t < theta_scan_min + 2*pi
            theta_t = theta_t - 2*pi;
        end
        % Map angle into the scan interval

        if theta_t < theta_scan_min || theta_t > theta_scan_max
            continue;
        end
        % Skip targets outside the scan range

        [~, idx_angle] = min(abs(theta_scan - theta_t));
        % Nearest scan angle index

        angle_window_deg = 10;
        % Angular spreading window

        angle_window_idx = round(angle_window_deg / rad2deg(theta_scan(2)-theta_scan(1)));
        % Convert angle window to index width

        for i = max(1, idx_angle-angle_window_idx):min(Nd, idx_angle+angle_window_idx)
            Pr_polar(r_idx(t), i) = Pr_polar(r_idx(t), i) + Pr_each(t, i);
        end
        % Accumulate power near the true target direction
    end

    Pr_polar_dBW = 10*log10(Pr_polar + realmin);
    % Convert to dBW

    [ThetaGrid, RGrid] = meshgrid(theta_scan, r_centers);
    % Build polar grid

    Xgrid = xr + RGrid .* cos(ThetaGrid);
    % Convert to x coordinates

    Ygrid = yr + RGrid .* sin(ThetaGrid);
    % Convert to y coordinates

    mask = (Xgrid >= x_min) & (Xgrid <= x_max) & ...
           (Ygrid >= y_min) & (Ygrid <= y_max);
    % Keep points inside the display window

    Xplot  = Xgrid(mask);
    Yplot  = Ygrid(mask);
    Prplot = Pr_polar_dBW(mask);
    % Output plot coordinates and power values

    maxVal = max(Prplot(:));
    % Maximum plotted value
end

function draw_scene(RectF1, RectR1, BikeR, ...
                    XF1, YF1, XR1, YR1, XR2, YR2, ...
                    XshF, YshF, XshR, YshR, ...
                    ant_xF, ant_yF, ant_xR, ant_yR, ...
                    Erect, xrF, yrF, xrR, yrR)
% Draw the scene geometry only.

    plot(RectF1.x, RectF1.y, 'k-', 'LineWidth', 2.0);
    % Draw front vehicle outline

    plot(RectR1.x, RectR1.y, 'k-', 'LineWidth', 2.0);
    % Draw rear vehicle outline

    plot(BikeR.x, BikeR.y, 'b-', 'LineWidth', 2.4);
    % Draw right-side bicycle line target

    plot(XF1, YF1, 'ko', 'MarkerFaceColor','y', 'MarkerSize',5);
    % Draw front vehicle scattering points

    plot(XR1, YR1, 'ko', 'MarkerFaceColor','c', 'MarkerSize',5);
    % Draw rear vehicle scattering points

    plot(XR2, YR2, 'bd', 'MarkerFaceColor','c', 'MarkerSize',5);
    % Draw bicycle scattering points

    % Front shadow points
    plot(XshF, YshF, 'rx', 'MarkerSize', 8, 'LineWidth', 1.2);
    % Mark front shadowed points

    % Rear shadow points allowed for display
    plot(XshR, YshR, 'rx', 'MarkerSize', 8, 'LineWidth', 1.2);
    % Mark rear shadowed points

    plot(ant_xF, ant_yF, 'm.', 'MarkerSize', 6);
    % Draw front array elements

    plot(ant_xR, ant_yR, 'm.', 'MarkerSize', 6);
    % Draw rear array elements

    plot(Erect.x, Erect.y, 'm-', 'LineWidth', 2.5);
    % Draw ego vehicle outline

    plot(xrF, yrF, 'mp', 'MarkerFaceColor','m', 'MarkerSize',12);
    % Draw front radar phase centre

    plot(xrR, yrR, 'mp', 'MarkerFaceColor','m', 'MarkerSize',12);
    % Draw rear radar phase centre
end