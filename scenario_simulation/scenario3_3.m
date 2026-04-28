clc; clear; close all;
% clc        : Clear previous text in the command window
% clear      : Clear existing variables in the workspace to avoid interference
% close all  : Close all open figure windows so the new plots are fresh

%% =========================
% Radar / array parameters
% =========================
c = 3e8;
% Speed of light in m/s

f0 = 25e9;
% Centre operating frequency of the radar / antenna system, in Hz

lambda0 = c / f0;
% Wavelength, lambda = c / f

d = lambda0/2;
% Spacing between adjacent array elements
% Here half-wavelength spacing is used

N = 20;
% Number of antenna elements

Pt_total = 100;
% Total transmit power of the array

Pt_elem  = Pt_total / N;
% Transmit power assigned to each element

n = 0:N-1;
% Element index, from 0 to N-1

w = 0.54 - 0.46*cos(2*pi*n/(N-1));
% Hamming taper
% Used to reduce sidelobes

w = w / sqrt(mean(w.^2));
% Normalize the taper weights

mC = (0:N-1) - (N-1)/2;
% Centred element index
% Used to place the array around its phase centre

%% =========================
% Ego vehicle
% =========================
L0 = 4.762;
% Ego vehicle length, in metres

W0 = 1.847;
% Ego vehicle width, in metres

xc0 = 0; yc0 = 0;
% Ego vehicle centre position

xL0 = xc0 - W0/2;
xR0 = xc0 + W0/2;
yB0 = yc0 - L0/2;
yT0 = yc0 + L0/2;
% Ego vehicle left, right, bottom, and top boundaries

Erect.x = [xL0 xR0 xR0 xL0 xL0];
Erect.y = [yB0 yB0 yT0 yT0 yB0];
% Ego vehicle rectangular outline

% Four radar phase centres
xrF  = xc0;  yrF  = yT0;
xrR  = xc0;  yrR  = yB0;
xrL  = xL0;  yrL  = yc0;
xrRt = xR0;  yrRt = yc0;
% Phase centres of the front, rear, left, and right arrays

% Element coordinates of the four arrays
ant_xF  = xrF  + mC*d;               ant_yF  = yrF  + zeros(size(mC));
ant_xR  = xrR  + mC*d;               ant_yR  = yrR  + zeros(size(mC));
ant_xL  = xrL  + zeros(size(mC));    ant_yL  = yrL  + mC*d;
ant_xRt = xrRt + zeros(size(mC));    ant_yRt = yrRt + mC*d;
% Front and rear arrays are arranged along x
% Left and right arrays are arranged along y

%% =========================
% Element pattern
% =========================
% Front / Rear use CIR pattern
T_FR = readtable('CIR_Gain_Phi90.xlsx');
% Read the antenna pattern data for the front / rear arrays

theta_all_FR = T_FR{:,1};
phi_all_FR   = T_FR{:,2};
gain_dBi_FR  = T_FR{:,3};
% Column 1: theta
% Column 2: phi
% Column 3: gain in dBi

mask_FR = abs(phi_all_FR - 90) < 1e-6 & theta_all_FR >= 0 & theta_all_FR <= 180;
% Extract the phi = 90 degree cut with theta from 0 to 180 degrees

theta_deg_FR = theta_all_FR(mask_FR);
G_lin_FR     = 10.^(gain_dBi_FR(mask_FR)/10);
% Convert gain from dBi to linear scale

[theta_deg_FR, idx_FR] = sort(theta_deg_FR);
G_lin_FR = G_lin_FR(idx_FR);
% Sort angles for stable interpolation

[theta_deg_FR, idxu_FR] = unique(theta_deg_FR, 'stable');
G_lin_FR = G_lin_FR(idxu_FR);
% Remove repeated angle values

element_gain_FR = @(a_rad) interp1(theta_deg_FR, G_lin_FR, abs(rad2deg(a_rad)), 'linear', 0);
% Element gain interpolation function for the front / rear arrays

% Left / Right use SQU pattern
T_LR = readtable('SQU_Gain_Theta25.xlsx');
% Read the antenna pattern data for the left / right arrays

theta_all_LR = T_LR{:,1};
phi_all_LR   = T_LR{:,2};
gain_dBi_LR  = T_LR{:,3};

mask_LR = abs(theta_all_LR - 25) < 1e-6 & phi_all_LR >= 0 & phi_all_LR <= 180;
% Extract the theta = 25 degree cut with phi from 0 to 180 degrees

phi_deg_LR = phi_all_LR(mask_LR);
G_lin_LR   = 10.^(gain_dBi_LR(mask_LR)/10);
% Convert gain from dBi to linear scale

[phi_deg_LR, idx_LR] = sort(phi_deg_LR);
G_lin_LR = G_lin_LR(idx_LR);
% Sort angles for stable interpolation

[phi_deg_LR, idxu_LR] = unique(phi_deg_LR, 'stable');
G_lin_LR = G_lin_LR(idxu_LR);
% Remove repeated angle values

element_gain_LR = @(a_rad) interp1(phi_deg_LR, G_lin_LR, abs(rad2deg(a_rad)), 'linear', 0);
% Element gain interpolation function for the left / right arrays

%% =========================
% Multi-frequency averaging
% =========================
Nf = 11;
% Number of frequency samples

B  = 400e6;
% System bandwidth, 400 MHz

f_vec      = f0 + linspace(-B/2, B/2, Nf);
% Frequency samples across the bandwidth

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
% Bicycle length in the simplified model

step = 0.2;
% Sampling interval along target edges

%% =========================
% Targets
% =========================
% Target 1: rear car
[XR1, YR1, RectR1, InfoR1] = make_rect_edges(-6.5, -12.0, Lcar, Wcar, step, 'RT');
% Rear vehicle
% 'RT' means right edge + top edge are used as scattering edges

% Target 2: front bicycle
[XF1, YF1, BikeF] = make_bicycle_line(0.0, 12.0, Lbik, step, 'vertical');
% Front bicycle, simplified as a vertical line target

% Target 3: right bicycle
[XRt1, YRt1, BikeRt] = make_bicycle_line(6.8, -1.5, Lbik, step, 'vertical');
% Right-side bicycle

% Target 4: right-front overlap car
[XRF1, YRF1, RectRF1, InfoRF1] = make_rect_edges(4.2, 4.0, Lcar, Wcar, step, 'LB');
% Right-front overlap vehicle
% 'LB' means left edge + bottom edge are used as scattering edges

% Combine all target scattering points
XptsAll = [XR1 XF1 XRt1 XRF1];
YptsAll = [YR1 YF1 YRt1 YRF1];

% 1 = rear car
% 2 = front bicycle
% 3 = right bicycle
% 4 = overlap car
typeAll = [ones(1,numel(XR1)), ...
           2*ones(1,numel(XF1)), ...
           3*ones(1,numel(XRt1)), ...
           4*ones(1,numel(XRF1))];

sigmaAll = ones(size(XptsAll));
% Set the RCS of each scattering point to 1

rng(1);
% Fix the random seed for repeatable results

phi_rand_all = 2*pi*rand(size(XptsAll));
% Random initial phase for each scattering point

%% =========================
% Scan angles
% =========================
Nd = 1081;
% Number of scan angle samples

theta_front = linspace(deg2rad(15),  deg2rad(165), Nd);
% Front scan range: 15 to 165 degrees

theta_rear  = linspace(deg2rad(195), deg2rad(345), Nd);
% Rear scan range: 195 to 345 degrees

theta_left  = linspace(deg2rad(150), deg2rad(210), Nd);
% Left scan range: 150 to 210 degrees

theta_right = linspace(deg2rad(-32), deg2rad(32), Nd);
% Right scan range: -32 to 32 degrees

%% =========================
% Visibility / shadow masks
% =========================
tol_vis = 1e-12;
% Numerical tolerance

visibleF  = (YptsAll >= yrF  - tol_vis);
visibleR  = (YptsAll <= yrR  + tol_vis);
visibleL  = (XptsAll <= xrL  + tol_vis);
visibleRt = (XptsAll >= xrRt - tol_vis);
% Check whether each point lies in the visible half-space of each array

shadowF  = ~visibleF;
shadowR  = ~visibleR;
shadowL  = ~visibleL;
shadowRt = ~visibleRt;
% Points outside the visible half-space are treated as shadowed

fprintf('\n');
fprintf('================ Visible Points by Array ================\n');
fprintf('Front array : %d / %d\n', sum(visibleF),  numel(visibleF));
fprintf('Rear  array : %d / %d\n', sum(visibleR),  numel(visibleR));
fprintf('Left  array : %d / %d\n', sum(visibleL),  numel(visibleL));
fprintf('Right array : %d / %d\n', sum(visibleRt), numel(visibleRt));
fprintf('=========================================================\n\n');

%% =========================
% Echo power
% =========================
Pr_each_F = compute_Pr_each_FR( ...
    XptsAll, YptsAll, sigmaAll, shadowF, xrF, yrF, ant_xF, ant_yF, ...
    theta_front, d, N, Pt_total, Pt_elem, w, phi_rand_all, ...
    k_vec, lambda_vec, element_gain_FR, +1);
% Received power from the front array
% +1 means the front array faces +y

Pr_each_R = compute_Pr_each_FR( ...
    XptsAll, YptsAll, sigmaAll, shadowR, xrR, yrR, ant_xR, ant_yR, ...
    theta_rear, d, N, Pt_total, Pt_elem, w, phi_rand_all, ...
    k_vec, lambda_vec, element_gain_FR, -1);
% Received power from the rear array
% -1 means the rear array faces -y

Pr_each_L = compute_Pr_each_LR( ...
    XptsAll, YptsAll, sigmaAll, shadowL, xrL, yrL, ant_xL, ant_yL, ...
    theta_left, d, N, Pt_total, Pt_elem, w, phi_rand_all, ...
    k_vec, lambda_vec, element_gain_LR, 'left');
% Received power from the left array

Pr_each_Rt = compute_Pr_each_LR( ...
    XptsAll, YptsAll, sigmaAll, shadowRt, xrRt, yrRt, ant_xRt, ant_yRt, ...
    theta_right, d, N, Pt_total, Pt_elem, w, phi_rand_all, ...
    k_vec, lambda_vec, element_gain_LR, 'right');
% Received power from the right array

%% =========================
% Target-wise contribution breakdown
% =========================
idxR1  = (typeAll == 1);
idxF1  = (typeAll == 2);
idxRt1 = (typeAll == 3);
idxRF1 = (typeAll == 4);
% Indices for the four target groups

print_target_breakdown('Rear Audi A4', idxR1,  Pr_each_F, Pr_each_R, Pr_each_L, Pr_each_Rt);
print_target_breakdown('Front bicycle', idxF1,  Pr_each_F, Pr_each_R, Pr_each_L, Pr_each_Rt);
print_target_breakdown('Right bicycle', idxRt1, Pr_each_F, Pr_each_R, Pr_each_L, Pr_each_Rt);
print_target_breakdown('Right-front overlap car', idxRF1, Pr_each_F, Pr_each_R, Pr_each_L, Pr_each_Rt);
% Print contribution breakdown for each target

%% =========================
% Polar accumulation
% =========================
r_step = 0.2;
% Range resolution of the polar grid

x_min = -15; x_max = 15;
y_min = -20; y_max = 20;
% Display range in the x-y plane

corners = [x_min y_min;
           x_min y_max;
           x_max y_min;
           x_max y_max];
% Corner coordinates of the display window

[XplotF, YplotF, PplotF, ~] = build_polar_map( ...
    XptsAll, YptsAll, Pr_each_F, xrF, yrF, theta_front, ...
    r_step, x_min, x_max, y_min, y_max, corners);
% Convert front-array polar map to x-y plot points

[XplotR, YplotR, PplotR, ~] = build_polar_map( ...
    XptsAll, YptsAll, Pr_each_R, xrR, yrR, theta_rear, ...
    r_step, x_min, x_max, y_min, y_max, corners);
% Convert rear-array polar map to x-y plot points

[XplotL, YplotL, PplotL, ~] = build_polar_map( ...
    XptsAll, YptsAll, Pr_each_L, xrL, yrL, theta_left, ...
    r_step, x_min, x_max, y_min, y_max, corners);
% Convert left-array polar map to x-y plot points

[XplotRt, YplotRt, PplotRt, ~] = build_polar_map( ...
    XptsAll, YptsAll, Pr_each_Rt, xrRt, yrRt, theta_right, ...
    r_step, x_min, x_max, y_min, y_max, corners);
% Convert right-array polar map to x-y plot points

marker_size = 28;
% Scatter marker size

global_max = max([PplotF(:); PplotR(:); PplotL(:); PplotRt(:)]);
% Maximum power value across all four arrays

dyn_range = 100;
% Display dynamic range for Figure 2

%% =========================
% Figure 1: front + rear only
% =========================
figure(1);
clf;
hold on;
axis equal;
grid on;
box on;

set(gcf, 'Color', 'w');
% Set figure background to white

cmap = jet(256);
colormap(cmap);

ax = gca;
ax.Color = 'none';
% Keep uncovered areas as the white background

plot([x_min x_max], [0 0], 'k-', 'LineWidth', 0.8);
plot([0 0], [y_min y_max], 'k-', 'LineWidth', 0.8);
% Draw x/y reference axes

xlim([x_min x_max]);
ylim([y_min y_max]);
% Set plot limits

scatter(XplotF, YplotF, marker_size, PplotF, 's', 'filled', 'MarkerEdgeColor', 'none');
scatter(XplotR, YplotR, marker_size, PplotR, 's', 'filled', 'MarkerEdgeColor', 'none');
% Plot front and rear array power distributions only

draw_scene(RectR1, BikeF, BikeRt, RectRF1, ...
           XR1, YR1, XF1, YF1, XRt1, YRt1, XRF1, YRF1, ...
           ant_xF, ant_yF, ant_xR, ant_yR, ant_xL, ant_yL, ant_xRt, ant_yRt, ...
           Erect, xrF, yrF, xrR, yrR, xrL, yrL, xrRt, yrRt, ...
           false, false);
% Overlay targets, scattering points, antenna elements, ego vehicle, and radar phase centres

caxis auto;
% Use automatic colour scaling for Figure 1

cb1 = colorbar;
ylabel(cb1, 'Range-angle cell P_r (dBW)');

title('Scenario 3 - Figure 1');
xlabel('x (m)');
ylabel('y (m)');
hold off;

%% =========================
% Figure 2: full fusion
% =========================
figure(2);
clf;
hold on;
axis equal;
grid on;
box on;

set(gcf, 'Color', 'w');
colormap(cmap);

ax = gca;
ax.Color = cmap(1,:);
% Set the axes background to the lowest colour in the colormap

plot([x_min x_max], [0 0], 'k-', 'LineWidth', 0.8);
plot([0 0], [y_min y_max], 'k-', 'LineWidth', 0.8);

xlim([x_min x_max]);
ylim([y_min y_max]);

scatter(XplotF,  YplotF,  marker_size, PplotF,  's', 'filled', 'MarkerEdgeColor', 'none');
scatter(XplotR,  YplotR,  marker_size, PplotR,  's', 'filled', 'MarkerEdgeColor', 'none');
scatter(XplotL,  YplotL,  marker_size, PplotL,  's', 'filled', 'MarkerEdgeColor', 'none');
scatter(XplotRt, YplotRt, marker_size, PplotRt, 's', 'filled', 'MarkerEdgeColor', 'none');
% Plot power distributions from all four arrays

draw_scene(RectR1, BikeF, BikeRt, RectRF1, ...
           XR1, YR1, XF1, YF1, XRt1, YRt1, XRF1, YRF1, ...
           ant_xF, ant_yF, ant_xR, ant_yR, ant_xL, ant_yL, ant_xRt, ant_yRt, ...
           Erect, xrF, yrF, xrR, yrR, xrL, yrL, xrRt, yrRt, ...
           true, true);
% Draw all four arrays in Figure 2

caxis([global_max-dyn_range, global_max]);
% Set the colour range from the maximum value down by 100 dB

cb2 = colorbar;
ylabel(cb2, 'Received power, P_r (dBW)');

title('Scenario 3');
xlabel('Lateral position, x (m)');
ylabel('Longitudinal position, y (m)');
hold off;

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
% rect   : full rectangle outline for plotting
% info   : rectangle boundary information

    build = @(a,b) linspace(a, b, round((b-a)/step)+1);
    % Generate points between a and b using approximately the given step

    xL = xc - W/2; xR = xc + W/2;
    % Left and right boundaries

    yB = yc - L/2; yT = yc + L/2;
    % Bottom and top boundaries

    switch upper(mode)
        case 'LB'
            xx = [xL*ones(size(build(yB,yT))), build(xL,xR)];
            yy = [build(yB,yT), yB*ones(size(build(xL,xR)))];
            % Left edge + bottom edge

        case 'RB'
            xx = [xR*ones(size(build(yB,yT))), build(xL,xR)];
            yy = [build(yB,yT), yB*ones(size(build(xL,xR)))];
            % Right edge + bottom edge

        case 'LT'
            xx = [xL*ones(size(build(yB,yT))), build(xL,xR)];
            yy = [build(yB,yT), yT*ones(size(build(xL,xR)))];
            % Left edge + top edge

        case 'RT'
            xx = [xR*ones(size(build(yB,yT))), build(xL,xR)];
            yy = [build(yB,yT), yT*ones(size(build(xL,xR)))];
            % Right edge + top edge

        otherwise
            error('Unknown mode');
            % Report an error for unknown edge modes
    end

    XY = unique([xx(:), yy(:)], 'rows', 'stable');
    % Combine coordinates and remove duplicated corner points

    X = XY(:,1).';
    Y = XY(:,2).';
    % Output as row vectors

    rect.x = [xL xR xR xL xL];
    rect.y = [yB yB yT yT yB];
    % Full rectangle outline for plotting

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
% xc, yc      : bicycle centre
% Lbik        : bicycle length
% step        : sampling interval
% orientation : horizontal / vertical / diag1 / diag2
%
% Outputs:
% X, Y        : sampled points on the bicycle line
% bikeLine    : line data for plotting

    s = -Lbik/2:step:Lbik/2;
    % Generate points along the bicycle length

    if abs(s(end)-Lbik/2) > 1e-9
        s(end+1) = Lbik/2;
        % Add the end point if it is not included
    end

    switch lower(orientation)
        case 'vertical'
            X = xc + zeros(size(s));
            Y = yc + s;
            % Vertical line: x is fixed, y changes

        case 'horizontal'
            X = xc + s;
            Y = yc + zeros(size(s));
            % Horizontal line: y is fixed, x changes

        case 'diag1'
            X = xc + s/sqrt(2);
            Y = yc + s/sqrt(2);
            % Diagonal line in the upper-right direction

        case 'diag2'
            X = xc + s/sqrt(2);
            Y = yc - s/sqrt(2);
            % Diagonal line in the lower-right direction

        otherwise
            error('Unknown bicycle orientation');
    end

    bikeLine.x = X;
    bikeLine.y = Y;
    % Store line data for plotting
end

function [Xplot, Yplot, Prplot, maxVal] = build_polar_map( ...
    Xpts, Ypts, Pr_each, xr, yr, theta_scan, ...
    r_step_cell, x_min, x_max, y_min, y_max, corners)
% Accumulate target-point power into a polar range-angle map,
% then map it back to x-y coordinates for plotting.

    Nd = numel(theta_scan);
    % Number of scan angle samples

    Rmax = max(hypot(corners(:,1)-xr, corners(:,2)-yr));
    % Maximum required display range from the radar

    Nr = floor(Rmax / r_step_cell);
    % Number of range cells

    r_centers = ((1:Nr) - 0.5) * r_step_cell;
    % Range-cell centre positions

    Rpts = hypot(Xpts - xr, Ypts - yr);
    % Distance from each target point to the radar

    r_idx = floor(Rpts / r_step_cell) + 1;
    % Range-cell index for each target point

    r_idx = max(1, min(Nr, r_idx));
    % Prevent index overflow

    theta_actual = atan2(Ypts - yr, Xpts - xr);
    % Actual angle of each target point relative to the radar

    theta_scan_min = min(theta_scan);
    theta_scan_max = max(theta_scan);
    % Scan angle limits

    Pr_polar = zeros(Nr, Nd);
    % Polar power map: rows are range cells, columns are angle cells

    for t = 1:numel(Xpts)
        theta_t = theta_actual(t);
        % Actual angle of the current target point

        if theta_t < theta_scan_min
            theta_t = theta_t + 2*pi;
        end
        if theta_t > theta_scan_max && theta_t < theta_scan_min + 2*pi
            theta_t = theta_t - 2*pi;
        end
        % Map angle into the same interval as the scan range

        if theta_t < theta_scan_min || theta_t > theta_scan_max
            continue;
        end
        % Skip points outside the scan range

        [~, idx_angle] = min(abs(theta_scan - theta_t));
        % Find the nearest scan angle index

        angle_window_deg = 10;
        % Angular spreading window for each target point

        angle_window_idx = round(angle_window_deg / rad2deg(theta_scan(2)-theta_scan(1)));
        % Convert angular window from degrees to index width

        for i = max(1, idx_angle-angle_window_idx):min(Nd, idx_angle+angle_window_idx)
            Pr_polar(r_idx(t), i) = Pr_polar(r_idx(t), i) + Pr_each(t, i);
        end
        % Add target power only around its true angular direction
    end

    Pr_polar_dBW = 10*log10(Pr_polar + realmin);
    % Convert polar power map to dBW

    [ThetaGrid, RGrid] = meshgrid(theta_scan, r_centers);
    % Construct polar grid

    Xgrid = xr + RGrid .* cos(ThetaGrid);
    Ygrid = yr + RGrid .* sin(ThetaGrid);
    % Convert polar coordinates to Cartesian coordinates

    mask = (Xgrid >= x_min) & (Xgrid <= x_max) & ...
           (Ygrid >= y_min) & (Ygrid <= y_max);
    % Keep only points inside the display window

    Xplot  = Xgrid(mask);
    Yplot  = Ygrid(mask);
    Prplot = Pr_polar_dBW(mask);
    % Output x-y plot points and power values

    maxVal = max(Prplot(:));
    % Maximum value in this plotted map
end

function Pr_each = compute_Pr_each_FR( ...
    Xpts, Ypts, sigma, is_shadow, xr, yr, ant_x, ant_y, ...
    theta_scan, d, N, Pt_total, Pt_elem, w, phi_rand, ...
    k_vec, lambda_vec, element_gain, facing_sign)
% Calculate received power for the front / rear arrays.
% These arrays are arranged along x, so they use the same phase form.
%
% Output:
% Pr_each(t,i) = received power of target point t at scan angle i

    Nd = numel(theta_scan);
    % Number of scan angles

    Nf = numel(k_vec);
    % Number of frequency samples

    Nt = numel(Xpts);
    % Number of target points

    mC = (0:N-1) - (N-1)/2;
    % Centred element index

    Rpts = hypot(Xpts - xr, Ypts - yr);
    % Distance from each target point to the radar phase centre

    Piso = Pt_elem ./ (4*pi*Rpts.^2);
    % Isotropic single-element reference power density

    Pr_each = zeros(Nt, Nd);
    % Initialise output matrix

    for fi = 1:Nf
        kf     = k_vec(fi);
        lambda = lambda_vec(fi);
        % Wavenumber and wavelength at the current frequency

        dphi = kf * d * cos(theta_scan);
        % Phase difference for arrays arranged along x

        for i = 1:Nd
            phi_feed = mC * dphi(i);
            % Feed phase for each element at this scan angle

            for t = 1:Nt
                if is_shadow(t)
                    continue;
                    % Skip target points that are shadowed for this array
                end

                dx = Xpts(t) - ant_x;
                dy = Ypts(t) - ant_y;
                rk = hypot(dx, dy);
                % Distance from the target point to each element

                a_k = atan2(dx, facing_sign * dy);
                % Local angle relative to the array main viewing direction

                G_elem = element_gain(a_k);
                % Element gain in this direction

                Pden_k = G_elem .* Pt_elem ./ (4*pi*rk.^2);
                % Power density produced by each element at the target point

                Vk = w .* sqrt(Pden_k);
                % Convert power density to field amplitude for coherent summation

                phi = kf*rk + phi_feed + phi_rand(t);
                % Total phase: propagation phase + feed phase + random initial phase

                E_sum = sum(Vk .* exp(1j*phi));
                % Coherent summation of complex fields

                P_array = abs(E_sum)^2;
                % Equivalent array power

                G = P_array / (N * Piso(t));
                % Equivalent array gain relative to N isotropic elements

                Pr_each(t,i) = Pr_each(t,i) + ...
                    Pt_total * G^2 * lambda^2 * sigma(t) / ((4*pi)^3 * Rpts(t)^4);
                % Monostatic radar equation
            end
        end
    end

    Pr_each = Pr_each / Nf;
    % Average over all frequency samples
end

function Pr_each = compute_Pr_each_LR( ...
    Xpts, Ypts, sigma, is_shadow, xr, yr, ant_x, ant_y, ...
    theta_scan, d, N, Pt_total, Pt_elem, w, phi_rand, ...
    k_vec, lambda_vec, element_gain, side_flag)
% Calculate received power for the left / right arrays.
% These arrays are arranged along y, so their phase form differs from the front / rear arrays.
%
% side_flag:
%   'left'  = left-facing array
%   'right' = right-facing array

    Nd = numel(theta_scan);
    % Number of scan angles

    Nf = numel(k_vec);
    % Number of frequency samples

    Nt = numel(Xpts);
    % Number of target points

    mC = (0:N-1) - (N-1)/2;
    % Centred element index

    Rpts = hypot(Xpts - xr, Ypts - yr);
    % Distance from each target point to the radar phase centre

    Piso = Pt_elem ./ (4*pi*Rpts.^2);
    % Isotropic single-element reference power density

    Pr_each = zeros(Nt, Nd);
    % Initialise output matrix

    for fi = 1:Nf
        kf       = k_vec(fi);
        lambda_f = lambda_vec(fi);
        % Wavenumber and wavelength at the current frequency

        dphi = kf * d * sin(theta_scan);
        % Phase difference for arrays arranged along y

        for i = 1:Nd
            phi_feed = mC * dphi(i);
            % Feed phase for each element at this scan angle

            for t = 1:Nt
                if is_shadow(t)
                    continue;
                    % Skip target points that are shadowed for this array
                end

                dx = Xpts(t) - ant_x;
                dy = Ypts(t) - ant_y;
                rk = hypot(dx, dy);
                % Distance from the target point to each element

                switch lower(side_flag)
                    case 'left'
                        a_k = atan2(dy, -dx);
                        % Main viewing direction of the left array is -x
                    case 'right'
                        a_k = atan2(dy, dx);
                        % Main viewing direction of the right array is +x
                    otherwise
                        error('Unknown side_flag');
                end

                G_elem = element_gain(a_k);
                % Element gain in this direction

                Pden_k = G_elem .* Pt_elem ./ (4*pi*rk.^2);
                % Power density produced by each element at the target point

                Vk = w .* sqrt(Pden_k);
                % Convert power density to field amplitude

                phi = kf*rk + phi_feed + phi_rand(t);
                % Total phase: propagation phase + feed phase + random initial phase

                E_sum = sum(Vk .* exp(1j*phi));
                % Coherent summation of complex fields

                P_array = abs(E_sum)^2;
                % Equivalent array power

                G = P_array / (N * Piso(t));
                % Equivalent array gain relative to N isotropic elements

                Pr_each(t,i) = Pr_each(t,i) + ...
                    Pt_total * G^2 * lambda_f^2 * sigma(t) / ((4*pi)^3 * Rpts(t)^4);
                % Monostatic radar equation
            end
        end
    end

    Pr_each = Pr_each / Nf;
    % Average over all frequency samples
end

function print_target_breakdown(target_name, idx_target, ...
    Pr_each_F, Pr_each_R, Pr_each_L, Pr_each_Rt)
% Print the contribution breakdown of one target across the four arrays.
% This includes total power sum and maximum power.

    sum_F  = sum(Pr_each_F(idx_target,:),  'all');
    sum_R  = sum(Pr_each_R(idx_target,:),  'all');
    sum_L  = sum(Pr_each_L(idx_target,:),  'all');
    sum_Rt = sum(Pr_each_Rt(idx_target,:), 'all');
    % Total power sum from each array

    max_F  = max(Pr_each_F(idx_target,:),  [], 'all');
    max_R  = max(Pr_each_R(idx_target,:),  [], 'all');
    max_L  = max(Pr_each_L(idx_target,:),  [], 'all');
    max_Rt = max(Pr_each_Rt(idx_target,:), [], 'all');
    % Maximum power from each array

    sum_total = sum_F + sum_R + sum_L + sum_Rt;
    max_total = max([max_F, max_R, max_L, max_Rt]);
    % Total sum and overall maximum

    fprintf('================ Contribution Breakdown ================\n');
    fprintf('\n%s:\n', target_name);
    fprintf('  Front : sum = %.6e, max = %.6e\n', sum_F,  max_F);
    fprintf('  Rear  : sum = %.6e, max = %.6e\n', sum_R,  max_R);
    fprintf('  Left  : sum = %.6e, max = %.6e\n', sum_L,  max_L);
    fprintf('  Right : sum = %.6e, max = %.6e\n', sum_Rt, max_Rt);
    fprintf('  Total : sum = %.6e, max = %.6e\n', sum_total, max_total);
    fprintf('========================================================\n\n');
end

function draw_scene(RectR1, BikeF, BikeRt, RectRF1, ...
                    XR1, YR1, XF1, YF1, XRt1, YRt1, XRF1, YRF1, ...
                    ant_xF, ant_yF, ant_xR, ant_yR, ant_xL, ant_yL, ant_xRt, ant_yRt, ...
                    Erect, xrF, yrF, xrR, yrR, xrL, yrL, xrRt, yrRt, ...
                    draw_side_arrays, draw_side_centres)
% Draw the scene only.
% This function is for visualisation and does not affect power calculation.

    plot(RectR1.x,  RectR1.y,  'k-', 'LineWidth', 1.8);
    plot(RectRF1.x, RectRF1.y, 'k-', 'LineWidth', 1.8);
    % Draw vehicle outlines

    plot(BikeF.x,  BikeF.y,  'b-', 'LineWidth', 2.2);
    plot(BikeRt.x, BikeRt.y, 'b-', 'LineWidth', 2.2);
    % Draw bicycle line targets

    plot(XR1,   YR1,   'ko', 'MarkerFaceColor','c', 'MarkerSize',5);
    plot(XRF1,  YRF1,  'ks', 'MarkerFaceColor','y', 'MarkerSize',5);
    plot(XF1,   YF1,   'bd', 'MarkerFaceColor','y', 'MarkerSize',5);
    plot(XRt1,  YRt1,  'bd', 'MarkerFaceColor','c', 'MarkerSize',5);
    % Draw scattering points

    plot(RectR1.x,  RectR1.y,  'k-', 'LineWidth', 2.2);
    plot(RectRF1.x, RectRF1.y, 'k-', 'LineWidth', 2.2);
    plot(BikeF.x,   BikeF.y,   'b-', 'LineWidth', 2.5);
    plot(BikeRt.x,  BikeRt.y,  'b-', 'LineWidth', 2.5);
    % Redraw target outlines with thicker lines

    plot(ant_xF, ant_yF, 'm.', 'MarkerSize', 6);
    plot(ant_xR, ant_yR, 'm.', 'MarkerSize', 6);
    % Draw front and rear array elements

    if draw_side_arrays
        plot(ant_xL,  ant_yL,  'm.', 'MarkerSize', 6);
        plot(ant_xRt, ant_yRt, 'm.', 'MarkerSize', 6);
        % Draw left and right array elements when required
    end

    plot(Erect.x, Erect.y, 'm-', 'LineWidth', 2.5);
    % Draw ego vehicle outline

    plot(xrF, yrF, 'mp', 'MarkerFaceColor','m', 'MarkerSize',12);
    plot(xrR, yrR, 'mp', 'MarkerFaceColor','m', 'MarkerSize',12);
    % Draw front and rear radar phase centres

    if draw_side_centres
        plot(xrL,  yrL,  'mp', 'MarkerFaceColor','m', 'MarkerSize',12);
        plot(xrRt, yrRt, 'mp', 'MarkerFaceColor','m', 'MarkerSize',12);
        % Draw left and right radar phase centres when required
    end
end