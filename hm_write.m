clc; clear; close all;

%% =========================
%% Load data
%% =========================
data1 = readtable('read2.txt', 'VariableNamingRule', 'preserve');
Vq1  = data1{:,2};
Vqb1 = data1{:,3};   % GREEN curve

data2 = readtable('write2.txt', 'VariableNamingRule', 'preserve');
Vq2  = data2{:,2};
Vqb2 = data2{:,3};   % RED curve

% Curve definitions
x1 = Vq1;   y1 = Vqb1;   % right curve (green)
x2 = Vqb2;   y2 = Vq2;   % left curve (red)

%% =========================
%% Build x(y) interpolators
%% =========================

% Green curve (right boundary)
[y1s, iy1] = sort(y1);
x1s = x1(iy1);
[y1s, ~, ic1] = unique(y1s);
x1s = accumarray(ic1, x1s, [], @mean);

% Red curve (left boundary)
[y2s, iy2] = sort(y2);
x2s = x2(iy2);
[y2s, ~, ic2] = unique(y2s);
x2s = accumarray(ic2, x2s, [], @mean);

% Interpolators
c_right = @(y) interp1(y1s, x1s, y, 'linear', NaN); % green
c_left  = @(y) interp1(y2s, x2s, y, 'linear', NaN); % red

%% =========================
%% SNM calculation
%% =========================
best_s = 0;
best_sq = [0 0 0 0];

y_min = max(min(y1s), min(y2s));
y_max = min(max(y1s), max(y2s));

for by = linspace(y_min, y_max, 500)

    xL = c_left(by);    % left boundary (red)
    xR = c_right(by);   % right boundary (green)

    if isnan(xL) || isnan(xR) || xR <= xL
        continue;
    end

    s = xR - xL;   % candidate square size

    % Check top edge
    y_top = by + s;
    if y_top > y_max
        continue;
    end

    xL_top = c_left(y_top);
    xR_top = c_right(y_top);

    if isnan(xL_top) || isnan(xR_top)
        continue;
    end

    % Ensure square fits at top also
    if (xR_top - xL_top) < s
        continue;
    end

    if s > best_s
        best_s = s;
        best_sq = [xL, by, s, s];
    end
end

SNM = best_s;

fprintf('SNM = %.4f V\n', SNM);

%% =========================
%% Plot
%% =========================
figure; set(gcf,'Color','k');
ax = axes;
ax.Color='k'; ax.XColor='w'; ax.YColor='w'; ax.GridColor='w';
hold on; grid on; axis equal;

plot(x1, y1, 'g', 'LineWidth', 2); % green curve
plot(x2, y2, 'r', 'LineWidth', 2); % red curve

% SNM square
rectangle('Position', best_sq, 'EdgeColor','y','LineWidth',2.5);

xlabel('Left node','Color','w');
ylabel('Right node','Color','w');
title(sprintf('SNM = %.4f V', SNM),'Color','w');

legend('Curve 1','Curve 2','SNM square',...
       'TextColor','w','Color','k');