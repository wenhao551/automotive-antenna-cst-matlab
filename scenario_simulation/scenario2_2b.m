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

Pt_elem = Pt_total / N;
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

xL0 = xc0 - W0/2;
xR0 = xc0 + W0/2;
yB0 = yc0 - L0/2;
yT0 = yc0 + L0/2;
% Ego vehicle boundaries

Erect.x = [xL0 xR0 xR0 xL0 xL0];
Erect.y = [yB0 yB0 yT0 yT0 yB0];
% Rectangular outline of the ego vehicle

EgoInfo.x_left  = xL0;
EgoInfo.x_right = xR0;
EgoInfo.y_bot   = yB0;
EgoInfo.y_top   = yT0;
% Store ego vehicle boundary information

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
% Read front / rear array element pattern data

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
% Remove duplicated angle values

element_gain_FR = @(a_rad) interp1(theta_deg_FR, G_lin_FR, abs(rad2deg(a_rad)), 'linear', 0);
% Element gain interpolation function for front / rear arrays

% Left / Right use SQU pattern
T_LR = readtable('SQU_Gain_Theta25.xlsx');
% Read left / right array element pattern data

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
% Remove duplicated angle values

element_gain_LR = @(a_rad) interp1(phi_deg_LR, G_lin_LR, abs(rad2deg(a_rad)), 'linear', 0);
% Element gain interpolation function for left / right arrays

%% =========================
% Multi-frequency averaging
% =========================
Nf = 11;
% Number of frequency samples

B = 400e6;
% System bandwidth, 400 MHz

f_vec = f0 + linspace(-B/2, B/2, Nf);
% Frequency samples around the centre frequency

lambda_vec = c ./ f_vec;
% Wavelength at each frequency sample

k_vec = 2*pi ./ lambda_vec;
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
% Targets
% =========================
% Target 1: left-upper car
[XF1, YF1, RectF1, InfoF1] = make_rect_edges(-4.5, 12.0, Lcar, Wcar, step, 'RB');
% 'RB' means right edge + bottom edge are used as scattering edges

% Target 2: left-lower car
[XR1, YR1, RectR1, InfoR1] = make_rect_edges(-4.5, -12.0, Lcar, Wcar, step, 'RT');
% 'RT' means right edge + top edge are used as scattering edges

% Target 3: right bicycle
[XRt1, YRt1, BikeRt] = make_bicycle_line(10.0, 0.0, Lbik, step, 'vertical');
% Right-side bicycle simplified as a vertical line target

% Combine all target scattering points
XptsAll = [XF1 XR1 XRt1];
YptsAll = [YF1 YR1 YRt1];

% 1 = left-upper car
% 2 = left-lower car
% 3 = right bicycle
typeAll = [ ...
    ones(1,numel(XF1)), ...
    2*ones(1,numel(XR1)), ...
    3*ones(1,numel(XRt1))];

sigmaAll = ones(size(XptsAll));
% Set the RCS of each scattering point to 1

%% =========================
% Shadow
% =========================
shadowF  = false(size(XptsAll));
shadowR  = false(size(XptsAll));
shadowL  = false(size(XptsAll));
shadowRt = false(size(XptsAll));
% Shadow masks for the front, rear, left, and right arrays

candF = (YptsAll >= yrF);
% Candidate points in the front half-space

shadowF(candF) = mark_shadow_by_blocker_multi( ...
    xrF, yrF, XptsAll(candF), YptsAll(candF), {InfoF1, EgoInfo});
% Check front-array shadowing using the left-upper car and ego vehicle

shadowF(typeAll == 1) = false;
% Prevent the left-upper car from shadowing itself

candR = (YptsAll <= yrR);
% Candidate points in the rear half-space

shadowR(candR) = mark_shadow_by_blocker_multi( ...
    xrR, yrR, XptsAll(candR), YptsAll(candR), {InfoR1, EgoInfo});
% Check rear-array shadowing using the left-lower car and ego vehicle

shadowR(typeAll == 2) = false;
% Prevent the left-lower car from shadowing itself

candL = (XptsAll <= xrL);
% Candidate points in the left half-space

shadowL(candL) = mark_shadow_by_blocker_multi( ...
    xrL, yrL, XptsAll(candL), YptsAll(candL), {EgoInfo});
% Check left-array shadowing using the ego vehicle

candRt = (XptsAll >= xrRt);
% Candidate points in the right half-space

shadowRt(candRt) = mark_shadow_by_blocker_multi( ...
    xrRt, yrRt, XptsAll(candRt), YptsAll(candRt), {EgoInfo});
% Check right-array shadowing using the ego vehicle

shadowRt(typeAll == 3) = false;
% Prevent the right-side bicycle from shadowing itself

fprintf('Front shadowed points : %d / %d\n', sum(shadowF),  numel(XptsAll));
fprintf('Rear  shadowed points : %d / %d\n', sum(shadowR),  numel(XptsAll));
fprintf('Left  shadowed points : %d / %d\n', sum(shadowL),  numel(XptsAll));
fprintf('Right shadowed points : %d / %d\n', sum(shadowRt), numel(XptsAll));
% Print the number of shadowed points for each array

%% =========================
% Random phase
% =========================
rng(1);
% Fix random seed for repeatable results

phi_rand_all = 2*pi*rand(size(XptsAll));
% Random initial phase for all scattering points

%% =========================
% Scan angles
% =========================
Nd = 1081;
% Number of scan angle samples

theta_front = linspace(0, pi, Nd);
% Front scan range: 0 to 180 degrees

theta_rear  = linspace(pi, 2*pi, Nd);
% Rear scan range: 180 to 360 degrees

theta_left  = linspace(pi/2, 3*pi/2, Nd);
% Left scan range: 90 to 270 degrees

theta_right = linspace(-pi/2, pi/2, Nd);
% Right scan range: -90 to +90 degrees

%% =========================
% Echo power
% =========================
Pr_each_F = compute_Pr_each_FR( ...
    XptsAll, YptsAll, sigmaAll, shadowF, xrF, yrF, ant_xF, ant_yF, ...
    theta_front, d, N, Pt_total, Pt_elem, w, phi_rand_all, ...
    k_vec, lambda_vec, element_gain_FR, +1);
% Received power from the front array

Pr_each_R = compute_Pr_each_FR( ...
    XptsAll, YptsAll, sigmaAll, shadowR, xrR, yrR, ant_xR, ant_yR, ...
    theta_rear, d, N, Pt_total, Pt_elem, w, phi_rand_all, ...
    k_vec, lambda_vec, element_gain_FR, -1);
% Received power from the rear array

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
% Corner coordinates of the display area

[XplotF, YplotF, PplotF_dBW, maxF] = build_polar_map( ...
    XptsAll, YptsAll, Pr_each_F, xrF, yrF, theta_front, ...
    r_step, x_min, x_max, y_min, y_max, corners);
% Convert front-array polar data to x-y plot points

[XplotR, YplotR, PplotR_dBW, maxR] = build_polar_map( ...
    XptsAll, YptsAll, Pr_each_R, xrR, yrR, theta_rear, ...
    r_step, x_min, x_max, y_min, y_max, corners);
% Convert rear-array polar data to x-y plot points

[XplotL, YplotL, PplotL_dBW, maxL] = build_polar_map( ...
    XptsAll, YptsAll, Pr_each_L, xrL, yrL, theta_left, ...
    r_step, x_min, x_max, y_min, y_max, corners);
% Convert left-array polar data to x-y plot points

[XplotRt, YplotRt, PplotRt_dBW, maxRt] = build_polar_map( ...
    XptsAll, YptsAll, Pr_each_Rt, xrRt, yrRt, theta_right, ...
    r_step, x_min, x_max, y_min, y_max, corners);
% Convert right-array polar data to x-y plot points

%% =========================
% Suppress weak arc responses
% =========================
arc_cut_dB = 12;
% Weak response suppression threshold in dB

PplotF_dBW(PplotF_dBW < maxF - arc_cut_dB) = -Inf;
PplotR_dBW(PplotR_dBW < maxR - arc_cut_dB) = -Inf;
PplotL_dBW(PplotL_dBW < maxL - arc_cut_dB) = -Inf;
PplotRt_dBW(PplotRt_dBW < maxRt - arc_cut_dB) = -Inf;
% Suppress weak responses for each array map

%% =========================
% Map to x-y (separate, linear power)
% =========================
PplotF_lin  = zeros(size(PplotF_dBW));
PplotR_lin  = zeros(size(PplotR_dBW));
PplotL_lin  = zeros(size(PplotL_dBW));
PplotRt_lin = zeros(size(PplotRt_dBW));
% Initialize linear power arrays

maskF_valid  = isfinite(PplotF_dBW);
maskR_valid  = isfinite(PplotR_dBW);
maskL_valid  = isfinite(PplotL_dBW);
maskRt_valid = isfinite(PplotRt_dBW);
% Valid points after weak-response suppression

PplotF_lin(maskF_valid)   = 10.^(PplotF_dBW(maskF_valid)/10);
PplotR_lin(maskR_valid)   = 10.^(PplotR_dBW(maskR_valid)/10);
PplotL_lin(maskL_valid)   = 10.^(PplotL_dBW(maskL_valid)/10);
PplotRt_lin(maskRt_valid) = 10.^(PplotRt_dBW(maskRt_valid)/10);
% Convert dBW back to linear power for fusion

%% =========================
% Fuse on a unified x-y grid
% =========================
xy_step = 0.2;
% Resolution of the unified x-y grid

xg = x_min:xy_step:x_max;
yg = y_min:xy_step:y_max;
% Unified x-y grid coordinates

Nx = numel(xg);
Ny = numel(yg);
% Number of grid points in x and y

Pxy_sum = zeros(Ny, Nx);
% Initialize fused power grid

Pxy_sum = accumulate_xy_grid_max(Pxy_sum, XplotF(maskF_valid),  YplotF(maskF_valid),  PplotF_lin(maskF_valid),  x_min, y_min, xy_step, Nx, Ny);
Pxy_sum = accumulate_xy_grid_max(Pxy_sum, XplotR(maskR_valid),  YplotR(maskR_valid),  PplotR_lin(maskR_valid),  x_min, y_min, xy_step, Nx, Ny);
Pxy_sum = accumulate_xy_grid_max(Pxy_sum, XplotL(maskL_valid),  YplotL(maskL_valid),  PplotL_lin(maskL_valid),  x_min, y_min, xy_step, Nx, Ny);
Pxy_sum = accumulate_xy_grid_max(Pxy_sum, XplotRt(maskRt_valid), YplotRt(maskRt_valid), PplotRt_lin(maskRt_valid), x_min, y_min, xy_step, Nx, Ny);
% Fuse four array maps on the unified grid using maximum power

Pxy_dBW = 10*log10(Pxy_sum + realmin);
% Convert fused linear power to dBW

global_max = max(Pxy_dBW(:));
% Maximum value in the fused map

dyn_range = 25;
% Final display dynamic range

Pxy_sum(Pxy_dBW < global_max - dyn_range) = 0;
% Remove points below the display threshold

Pxy_dBW = 10*log10(Pxy_sum + realmin);
% Recalculate dBW after thresholding

[XG, YG] = meshgrid(xg, yg);
% Build full x-y grid

mask_plot = Pxy_sum > 0;
% Keep only valid plot points

XplotAll = XG(mask_plot);
YplotAll = YG(mask_plot);
PplotAll = Pxy_dBW(mask_plot);
% Extract final plot coordinates and power values

%% =========================
% Plot
% =========================
figure;
clf;
hold on;
axis equal;
grid on;
box on;

set(gcf, 'Color', 'w');
% Set figure background to white

cmap = jet(256);
colormap(cmap);
% Set colormap

ax = gca;
ax.Color = cmap(1,:);
% Set axes background to the lowest colormap colour

plot([x_min x_max], [0 0], 'k-', 'LineWidth', 0.8);
plot([0 0], [y_min y_max], 'k-', 'LineWidth', 0.8);
% Draw x and y reference axes

xlim([x_min x_max]);
ylim([y_min y_max]);
% Set display limits

scatter(XplotAll, YplotAll, 28, PplotAll, 's', 'filled', 'MarkerEdgeColor', 'none');
% Plot the final fused power map

draw_scene( ...
    RectF1, RectR1, BikeRt, ...
    XF1, YF1, XR1, YR1, XRt1, YRt1, ...
    XptsAll(shadowF | shadowR | shadowL | shadowRt), ...
    YptsAll(shadowF | shadowR | shadowL | shadowRt), ...
    ant_xF, ant_yF, ant_xR, ant_yR, ant_xL, ant_yL, ant_xRt, ant_yRt, ...
    Erect, xrF, yrF, xrR, yrR, xrL, yrL, xrRt, yrRt);
% Overlay scene geometry, scattering points, shadowed points, arrays, and ego vehicle

caxis([global_max-dyn_range, global_max]);
% Set colour range

cb = colorbar;
ylabel(cb, 'Received power, P_r (dBW)');

title('Scenario 2b');
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
% rect   : complete rectangle outline for plotting
% info   : rectangle boundary information

    build = @(a,b) linspace(a, b, round((b-a)/step)+1);
    % Generate points between a and b

    xL = xc - W/2; xR = xc + W/2;
    % Left and right boundaries

    yB = yc - L/2; yT = yc + L/2;
    % Bottom and top boundaries

    switch upper(mode)
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
    end

    XY = unique([xx(:), yy(:)], 'rows', 'stable');
    % Remove duplicated corner points

    X = XY(:,1).';
    Y = XY(:,2).';
    % Output x and y as row vectors

    rect.x = [xL xR xR xL xL];
    rect.y = [yB yB yT yT yB];
    % Complete rectangle outline

    info.x_left  = xL;
    info.x_right = xR;
    info.y_bot   = yB;
    info.y_top   = yT;
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

function shadow = mark_shadow_by_blocker_multi(xr, yr, X, Y, info_list)
% Check shadowing from multiple rectangular blockers.

    shadow = false(size(X));
    % Initial shadow mask

    for m = 1:numel(info_list)
        shadow = shadow | mark_shadow_by_blocker(xr, yr, X, Y, info_list{m});
        % A point is shadowed if it is blocked by any blocker
    end
end

function shadow = mark_shadow_by_blocker(xr, yr, X, Y, info)
% Check whether target points are blocked by one rectangular blocker.

    tol = 1e-9;
    % Numerical tolerance

    p1 = [info.x_right, info.y_bot];
    p2 = [info.x_left,  info.y_top];
    % Two diagonal points of the blocker

    v = p2 - p1;
    % Diagonal direction vector

    radar_vec = [xr - p1(1), yr - p1(2)];
    % Vector from p1 to radar

    radar_side = v(1)*radar_vec(2) - v(2)*radar_vec(1);
    % Side of the diagonal where the radar is located

    th1 = atan2(p1(2)-yr, p1(1)-xr);
    th2 = atan2(p2(2)-yr, p2(1)-xr);
    % Angular limits of the blocker

    if abs(th1-th2) > pi
        if th1 < th2
            th1 = th1 + 2*pi;
        else
            th2 = th2 + 2*pi;
        end
    end
    % Handle angle wrapping

    th_min = min(th1, th2);
    th_max = max(th1, th2);
    % Shadow angular range

    shadow = false(size(X));
    % Initial output mask

    for k = 1:numel(X)
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

function Pr_each = compute_Pr_each_FR( ...
    Xpts, Ypts, sigma, is_shadow, xr, yr, ant_x, ant_y, ...
    theta_scan, d, N, Pt_total, Pt_elem, w, phi_rand, ...
    k_vec, lambda_vec, element_gain, facing_sign)
% Calculate received power for front / rear arrays.
% Front and rear arrays are arranged along x.

    Nd = numel(theta_scan);
    Nf = numel(k_vec);
    Nt = numel(Xpts);

    mC = (0:N-1) - (N-1)/2;
    Rpts = hypot(Xpts - xr, Ypts - yr);
    Piso = Pt_elem ./ (4*pi*Rpts.^2);

    Pr_each = zeros(Nt, Nd);

    for fi = 1:Nf
        kf = k_vec(fi);
        lambda = lambda_vec(fi);

        dphi = kf * d * cos(theta_scan);
        % Phase difference for arrays arranged along x

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
                % Local angle relative to front or rear facing direction

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
    % Average over all frequency samples
end

function Pr_each = compute_Pr_each_LR( ...
    Xpts, Ypts, sigma, is_shadow, xr, yr, ant_x, ant_y, ...
    theta_scan, d, N, Pt_total, Pt_elem, w, phi_rand, ...
    k_vec, lambda_vec, element_gain, side_flag)
% Calculate received power for left / right arrays.
% Left and right arrays are arranged along y.

    Nd = numel(theta_scan);
    Nf = numel(k_vec);
    Nt = numel(Xpts);

    mC = (0:N-1) - (N-1)/2;
    Rpts = hypot(Xpts - xr, Ypts - yr);
    Piso = Pt_elem ./ (4*pi*Rpts.^2);

    Pr_each = zeros(Nt, Nd);

    for fi = 1:Nf
        kf       = k_vec(fi);
        lambda_f = lambda_vec(fi);

        dphi = kf * d * sin(theta_scan);
        % Phase difference for arrays arranged along y

        for i = 1:Nd
            phi_feed = mC * dphi(i);

            for t = 1:Nt
                if is_shadow(t)
                    continue;
                end

                dx = Xpts(t) - ant_x;
                dy = Ypts(t) - ant_y;
                rk = hypot(dx, dy);

                switch lower(side_flag)
                    case 'left'
                        a_k = atan2(dy, -dx);
                        % Left array faces -x
                    case 'right'
                        a_k = atan2(dy, dx);
                        % Right array faces +x
                    otherwise
                        error('Unknown side_flag');
                end

                G_elem = element_gain(a_k);
                Pden_k = G_elem .* Pt_elem ./ (4*pi*rk.^2);
                Vk = w .* sqrt(Pden_k);

                phi = kf*rk + phi_feed + phi_rand(t);
                E_sum = sum(Vk .* exp(1j*phi));
                P_array = abs(E_sum)^2;
                G = P_array / (N * Piso(t));

                Pr_each(t,i) = Pr_each(t,i) + ...
                    Pt_total * G^2 * lambda_f^2 * sigma(t) / ((4*pi)^3 * Rpts(t)^4);
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
    Ygrid = yr + RGrid .* sin(ThetaGrid);
    % Convert polar coordinates to Cartesian coordinates

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

function Pxy_sum = accumulate_xy_grid_max(Pxy_sum, X, Y, P, x_min, y_min, xy_step, Nx, Ny)
% Map scattered x-y power points to a regular grid.
% If multiple values fall into the same grid cell, keep the maximum value.

    ix = floor((X - x_min)/xy_step) + 1;
    iy = floor((Y - y_min)/xy_step) + 1;
    % Convert coordinates to grid indices

    valid = ix >= 1 & ix <= Nx & iy >= 1 & iy <= Ny;
    % Check valid indices

    ix = ix(valid);
    iy = iy(valid);
    P  = P(valid);
    % Keep valid points only

    for k = 1:numel(P)
        Pxy_sum(iy(k), ix(k)) = max(Pxy_sum(iy(k), ix(k)), P(k));
        % Keep the strongest response in each grid cell
    end
end

function draw_scene( ...
    RectF1, RectR1, BikeRt, ...
    XF1, YF1, XR1, YR1, XRt1, YRt1, ...
    Xsh, Ysh, ant_xF, ant_yF, ant_xR, ant_yR, ant_xL, ant_yL, ant_xRt, ant_yRt, ...
    Erect, xrF, yrF, xrR, yrR, xrL, yrL, xrRt, yrRt)
% Draw scene geometry only.

    plot(RectF1.x, RectF1.y, 'k-', 'LineWidth', 2.0);
    % Draw left-upper car outline

    plot(RectR1.x, RectR1.y, 'k-', 'LineWidth', 2.0);
    % Draw left-lower car outline

    plot(BikeRt.x, BikeRt.y, 'b-', 'LineWidth', 2.4);
    % Draw right-side bicycle line target

    plot(XF1, YF1, 'ko', 'MarkerFaceColor','y', 'MarkerSize',5);
    % Draw left-upper car scattering points

    plot(XR1, YR1, 'ko', 'MarkerFaceColor','c', 'MarkerSize',5);
    % Draw left-lower car scattering points

    plot(XRt1, YRt1, 'bd', 'MarkerFaceColor','g', 'MarkerSize',5);
    % Draw right-side bicycle scattering points

    plot(Xsh, Ysh, 'rx', 'MarkerSize', 8, 'LineWidth', 1.2);
    % Mark shadowed points

    plot(ant_xF,  ant_yF,  'm.', 'MarkerSize', 6);
    plot(ant_xR,  ant_yR,  'm.', 'MarkerSize', 6);
    plot(ant_xL,  ant_yL,  'm.', 'MarkerSize', 6);
    plot(ant_xRt, ant_yRt, 'm.', 'MarkerSize', 6);
    % Draw array element positions

    plot(Erect.x, Erect.y, 'm-', 'LineWidth', 2.5);
    % Draw ego vehicle outline

    plot(xrF,  yrF,  'mp', 'MarkerFaceColor','m', 'MarkerSize',12);
    plot(xrR,  yrR,  'mp', 'MarkerFaceColor','m', 'MarkerSize',12);
    plot(xrL,  yrL,  'mp', 'MarkerFaceColor','m', 'MarkerSize',12);
    plot(xrRt, yrRt, 'mp', 'MarkerFaceColor','m', 'MarkerSize',12);
    % Draw radar phase centres
end