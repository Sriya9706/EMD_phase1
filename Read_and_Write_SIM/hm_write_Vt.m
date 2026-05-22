clc; clear; close all;

%% ── Parse all MC runs ───────────────────────────────────────────────────────
function [x, y] = parse_mc_runs(filename)
    fid = fopen(filename, 'r');
    lines = {}; while ~feof(fid); line = strtrim(fgetl(fid)); lines{end+1} = line; end
    fclose(fid);
    x = {}; y = {}; current_run = 0; temp_x = []; temp_y = [];
    for i = 1:length(lines)
        line = lines{i};
        if contains(line, 'Step Information')
            if current_run > 0 && ~isempty(temp_x)
                x{current_run} = temp_x; y{current_run} = temp_y;
            end
            token = regexp(line, 'Mc_runs=(\d+)', 'tokens');
            current_run = str2double(token{1}{1});
            temp_x = []; temp_y = [];
        else
            nums = sscanf(line, '%f %f');
            if length(nums) == 2
                temp_x = [temp_x; nums(1)]; temp_y = [temp_y; nums(2)];
            end
        end
    end
    if current_run > 0 && ~isempty(temp_x); x{current_run} = temp_x; y{current_run} = temp_y; end
end

[y1_all, x1_all] = parse_mc_runs('lastread.txt');   % green (read)
[x2_all, y2_all] = parse_mc_runs('lastwrite.txt');  % red   (write)

num_runs = length(x1_all);

%% ── Build common x grid and interpolate all runs onto it ───────────────────
N = 5000;

% Green envelope: common x range across all green runs
xg_lo = max(cellfun(@min, x1_all));
xg_hi = min(cellfun(@max, x1_all));
xg    = linspace(xg_lo, xg_hi, N)';

Y_green = zeros(N, num_runs);
for k = 1:num_runs
    [xs, ix] = sort(x1_all{k}); [xs, ux] = unique(xs);
    ys = y1_all{k}(ix); ys = ys(ux);
    Y_green(:,k) = interp1(xs, ys, xg, 'linear', NaN);
end
y_green_env = min(Y_green, [], 2);   % innermost = minimum y

% Red envelope: common x range across all red runs
xr_lo = max(cellfun(@min, x2_all));
xr_hi = min(cellfun(@max, x2_all));
xr    = linspace(xr_lo, xr_hi, N)';

Y_red = zeros(N, num_runs);
for k = 1:num_runs
    [xs, ix] = sort(x2_all{k}); [xs, ux] = unique(xs);
    ys = y2_all{k}(ix); ys = ys(ux);
    Y_red(:,k) = interp1(xs, ys, xr, 'linear', NaN);
end
y_red_env = max(Y_red, [], 2);       % innermost = maximum y

%% ── SNM on envelope curves (hm_write_new logic) ────────────────────────────
% Green envelope: x1=xg, y1=y_green_env  → right boundary
% Red envelope:   x2=xr, y2=y_red_env    → left boundary

[y1s,iy1]=sort(y_green_env); x1s=xg(iy1); [y1s,u1]=unique(y1s); x1s=x1s(u1);
[y2s,iy2]=sort(y_red_env);   x2s=xr(iy2); [y2s,u2]=unique(y2s); x2s=x2s(u2);

c_green = @(y) interp1(y1s, x1s, y, 'linear', NaN);
c_red   = @(y) interp1(y2s, x2s, y, 'linear', NaN);

y_lo = max(y1s(1), y2s(1));
y_hi = min(y1s(end), y2s(end));
yg_snm = linspace(y_lo, y_hi, N)';
xG = c_green(yg_snm);
xR = c_red(yg_snm);

best_s = 0; best_sq = [0 0 0 0];
for i = 1:N
    yb = yg_snm(i); s_lo = 0; s_hi = y_hi - yb;
    for iter = 1:80
        s_mid = (s_lo + s_hi) / 2;
        mask  = yg_snm >= yb & yg_snm <= yb + s_mid;
        if ~any(mask), s_hi = s_mid; continue; end
        avail = min(xG(mask)) - max(xR(mask));
        if avail >= s_mid; s_lo = s_mid; else; s_hi = s_mid; end
    end
    s = s_lo;
    if s > best_s
        mask = yg_snm >= yb & yg_snm <= yb + s;
        best_s = s; best_sq = [max(xR(mask)), yb, s, s];
    end
end

fprintf('Worst-case SNM = %.4f V\n', best_s);

%% ── Plot ────────────────────────────────────────────────────────────────────
figure; set(gcf,'Color','k');
ax=axes; ax.Color='k'; ax.XColor='w'; ax.YColor='w'; ax.GridColor='w';
hold on; grid on; axis equal;

% All runs
for k = 1:num_runs
    plot(x1_all{k}, y1_all{k}, 'g', 'LineWidth', 0.8);
    plot(x2_all{k}, y2_all{k}, 'r', 'LineWidth', 0.8);
end


rectangle('Position', best_sq, 'EdgeColor', 'y', 'LineWidth', 2.5);

xlabel('Left node',  'Color', 'w');
ylabel('Right node', 'Color', 'w');
title(sprintf('Worst-case SNM = %.4f V', best_s), 'Color', 'w');
legend('Read runs','Write runs','SNM square','TextColor','w','Color','k');