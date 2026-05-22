clc; clear; close all;

%% =========================================================
%% Worst-Case SNM from Monte Carlo Butterfly Curves
%%
%% From the plot:
%%   Green = Curve 1: x=Vq,  y=Vqb  (raw data)
%%   Red   = Curve 2: x=Vqb, y=Vq   (same data, axes flipped)
%%
%% Lobe 1 (LEFT side of metastable point):
%%   Worst-case = tightest gap between green and red
%%   -> innermost green = pointwise MIN of all green curves (lowest green)
%%   -> innermost red   = pointwise MAX of all red curves   (highest red)
%%
%% Lobe 2 (RIGHT side of metastable point):
%%   -> innermost green = pointwise MAX of all green curves (highest green)
%%   -> innermost red   = pointwise MIN of all red curves   (lowest red)
%%
%% Then feed these two envelope curves into the exact same SNM
%% algorithm from Nosie_Margin_calculation_for_read_and_hold.m
%% =========================================================

%% ----------------------------------------------------------
%% 1.  Parse all MC runs
%% ----------------------------------------------------------
[Vq_all, Vqb_all] = parse_mc_runs('lastread.txt');
num_runs = length(Vq_all);
fprintf('Total MC runs: %d\n', num_runs);

%% ----------------------------------------------------------
%% 2.  Common interpolation grid
%% ----------------------------------------------------------
all_x = cell2mat(Vq_all(:));   % Vq and Vqb share the same range
x_min = min(all_x);
x_max = max(all_x);
N     = 3000;
xg    = linspace(x_min, x_max, N)';

%% ----------------------------------------------------------
%% 3.  Interpolate every run onto the common grid
%%     Green (C1): y = Vqb as function of x = Vq
%%     Red   (C2): y = Vq  as function of x = Vqb
%% ----------------------------------------------------------
G = nan(N, num_runs);   % green curves on grid
R = nan(N, num_runs);   % red   curves on grid

for k = 1:num_runs
    % Green: x=Vq, y=Vqb
    [xs, ix] = sort(Vq_all{k});
    ys = Vqb_all{k}(ix);
    [xs, ux] = unique(xs, 'stable'); ys = ys(ux);
    G(:,k) = interp1(xs, ys, xg, 'linear', NaN);

    % Red: x=Vqb, y=Vq
    [xs, ix] = sort(Vqb_all{k});
    ys = Vq_all{k}(ix);
    [xs, ux] = unique(xs, 'stable'); ys = ys(ux);
    R(:,k) = interp1(xs, ys, xg, 'linear', NaN);
end

%% ----------------------------------------------------------
%% 4.  Find global metastable point P2x
%% ----------------------------------------------------------
P2x_all = nan(1, num_runs);
for k = 1:num_runs
    [xi, yi] = polyxpoly(Vq_all{k}, Vqb_all{k}, Vqb_all{k}, Vq_all{k});
    if ~isempty(xi)
        [~, mid] = min(abs(xi - yi));
        P2x_all(k) = xi(mid);
    end
end
P2x = median(P2x_all, 'omitnan');
fprintf('Global metastable P2x = %.4f V\n', P2x);

[~, p2_idx] = min(abs(xg - P2x));   % grid index of P2x

%% ----------------------------------------------------------
%% 5.  Build innermost (worst-case) envelope curves
%%
%%  Lobe 1 (x <= P2x):
%%    innermost green = MIN across runs (lowest green squeezes lobe)
%%    innermost red   = MAX across runs (highest red squeezes lobe)
%%
%%  Lobe 2 (x >= P2x):
%%    innermost green = MAX across runs (highest green squeezes lobe)
%%    innermost red   = MIN across runs (lowest red squeezes lobe)
%% ----------------------------------------------------------
inner_green_L1 = min(G(1:p2_idx, :), [], 2, 'omitnan');  % left  of P2x
inner_red_L1   = max(R(1:p2_idx, :), [], 2, 'omitnan');

inner_green_L2 = max(G(p2_idx:end, :), [], 2, 'omitnan'); % right of P2x
inner_red_L2   = min(R(p2_idx:end, :), [], 2, 'omitnan');

xg_L1 = xg(1:p2_idx);
xg_L2 = xg(p2_idx:end);

%% ----------------------------------------------------------
%% 6.  Stitch into full-length x1,y1 / x2,y2 for the SNM code
%%     Curve 1 (green): use L1-inner on left,  L2-inner on right
%%     Curve 2 (red):   use L1-inner on left,  L2-inner on right
%% ----------------------------------------------------------
x1 = [xg_L1;          xg_L2];
y1 = [inner_green_L1; inner_green_L2];   % worst-case green (Vq vs Vqb)

x2 = [xg_L1;        xg_L2];
y2 = [inner_red_L1; inner_red_L2];       % worst-case red   (Vqb vs Vq)

%% Remove NaNs
valid = ~isnan(y1) & ~isnan(y2);
x1 = x1(valid); y1 = y1(valid);
x2 = x2(valid); y2 = y2(valid);

%% ----------------------------------------------------------
%% 7.  Run EXACT SNM algorithm from
%%     Nosie_Margin_calculation_for_read_and_hold.m
%% ----------------------------------------------------------

%% Find metastable point on the stitched curves
[xi, yi] = polyxpoly(x1, y1, x2, y2);
[~, mid] = min(abs(xi - yi));
P2x_snm = xi(mid);
P2y_snm = yi(mid);
fprintf('SNM metastable: (%.4f, %.4f)\n', P2x_snm, P2y_snm);

%% Interpolators
[x1s, ix1] = sort(x1); y1s = y1(ix1);
[x2s, ix2] = sort(x2); y2s = y2(ix2);
[x1s, ux1] = unique(x1s,'stable'); y1s = y1s(ux1);
[x2s, ux2] = unique(x2s,'stable'); y2s = y2s(ux2);

c1_y = @(x) interp1(x1s, y1s, x, 'linear', NaN);
c2_y = @(x) interp1(x2s, y2s, x, 'linear', NaN);


%% --- LOBE 1 (upper-left lobe: force into high-Vqb region) ---
best_s1  = 0;
best_sq1 = [0 0 0 0];
tol = 1e-3;

% From the plot: Lobe 1 is clearly in x ~ [0.05, 0.25], y ~ [0.6, 0.9]
% Search bx only in the LEFT portion, well away from metastable point
x_lobe1_max = P2x_snm * 0.5;   % only left half of left side
y_lobe1_min = P2y_snm + 0.15;  % must be well above metastable y

for bx = linspace(0, x_lobe1_max, 400)
    y1_at_bx = c1_y(bx);
    y2_at_bx = c2_y(bx);
    if isnan(y1_at_bx) || isnan(y2_at_bx), continue; end

    by = min(y1_at_bx, y2_at_bx);

    % Force into upper lobe
    if by < y_lobe1_min, continue; end

    for s = linspace(0.001, 0.5, 400)
        if bx + s > x_lobe1_max + 0.05, break; end

        y1_at_bxs = c1_y(bx + s);
        y2_at_bxs = c2_y(bx + s);
        if isnan(y1_at_bxs) || isnan(y2_at_bxs), continue; end

        top = by + s;

        % Top must not exceed either curve at left and right edges
        if top > max(y1_at_bx,  y2_at_bx)  + tol, continue; end
        if top > max(y1_at_bxs, y2_at_bxs) + tol, continue; end

        % Bottom-right must still be above lower curve
        if by < min(y1_at_bxs, y2_at_bxs) - tol, continue; end

        % All 4 corners must lie between the two curves
        corners_x = [bx, bx+s, bx,   bx+s];
        corners_y = [by, by,   by+s, by+s];
        valid_sq = true;
        for ci = 1:4
            lo = min(c1_y(corners_x(ci)), c2_y(corners_x(ci)));
            hi = max(c1_y(corners_x(ci)), c2_y(corners_x(ci)));
            if isnan(lo) || isnan(hi), valid_sq = false; break; end
            if corners_y(ci) > hi + tol || corners_y(ci) < lo - tol
                valid_sq = false; break;
            end
        end

        if valid_sq && s > best_s1
            best_s1  = s;
            best_sq1 = [bx, by, s, s];
        end
    end
end

%% --- LOBE 2 (right side) ---
best_s2  = 0;
best_sq2 = [0 0 0 0];

idx1 = find(x1 > P2x_snm);
idx2 = find(x2 > P2x_snm);

for i = idx1'
    x_bot = x1(i); y_bot = y1(i);
    for j = idx2'
        x_top = x2(j); y_top = y2(j);
        s = min(y_top - y_bot, x_top - x_bot);
        if s <= 0, continue; end
        xs = [x_bot, x_bot+s, x_bot, x_bot+s];
        ys = [y_bot, y_bot,   y_bot+s, y_bot+s];
        valid_sq = true;
        for k = 1:4
            xk = xs(k); yk = ys(k);
            y1k = c1_y(xk); y2k = c2_y(xk);
            if isnan(y1k) || isnan(y2k), valid_sq=false; break; end
            if yk > max(y1k,y2k)+tol || yk < min(y1k,y2k)-tol
                valid_sq=false; break;
            end
        end
        if valid_sq && s > best_s2
            best_s2  = s;
            best_sq2 = [x_bot, y_bot, s, s];
        end
    end
end

%% ----------------------------------------------------------
%% 8.  Final SNM
%% ----------------------------------------------------------
SNM = min(best_s1, best_s2);

fprintf('\nLobe 1 SNM = %.4f V\n', best_s1);
fprintf('Lobe 2 SNM = %.4f V\n', best_s2);
fprintf('Worst-Case SNM = %.4f V\n', SNM);

%% ----------------------------------------------------------
%% 9.  Plot
%% ----------------------------------------------------------
figure; set(gcf,'Color','k','Position',[100 100 900 700]);
ax = axes;
ax.Color='k'; ax.XColor='w'; ax.YColor='w'; ax.GridColor='w';
hold on; grid on; axis equal;

%% All MC runs — same intensity, no transparency
for k = 1:num_runs
    plot(Vq_all{k},  Vqb_all{k}, 'Color',[0.1 0.6 0.1], 'LineWidth', 0.8);
    plot(Vqb_all{k}, Vq_all{k},  'Color',[0.6 0.1 0.1], 'LineWidth', 0.8);
end

%% Worst-case envelope — same color, slightly thicker to stay visible
plot(x1, y1, 'Color',[0.1 0.6 0.1], 'LineWidth', 1.5);
plot(x2, y2, 'Color',[0.6 0.1 0.1], 'LineWidth', 1.5);

%% y = x diagonal
lim = x_max * 1.05;
plot([0 lim],[0 lim], 'w--', 'LineWidth', 1.5);

%% Squares
rectangle('Position', best_sq1, 'EdgeColor','y', 'LineWidth', 2.5);
rectangle('Position', best_sq2, 'EdgeColor','c', 'LineWidth', 2.5);

%% Metastable point
plot(xi, yi, 'ro', 'MarkerSize', 8, 'LineWidth', 2);

xlabel('Left node (Vq)',   'Color','w');
ylabel('Right node (Vqb)', 'Color','w');
title(sprintf('Worst-Case SNM = %.4f V   (Lobe1=%.4f V,  Lobe2=%.4f V)', ...
    SNM, best_s1, best_s2), 'Color','w');
legend('MC green','MC red','Innermost green','Innermost red','y=x', ...
       'Lobe1 sq','Lobe2 sq','Intersections', ...
       'TextColor','w','Color','k','Location','northeast');