clc; clear; close all;

%% Load data
data = readtable('6T_cell.txt', 'VariableNamingRule', 'preserve');
Vq  = data{:,2};
Vqb = data{:,3};

x1 = Vq;   y1 = Vqb;   % Curve 1
x2 = Vqb;  y2 = Vq;    % Curve 2

%% Find metastable point (middle intersection)
[xi, yi] = polyxpoly(x1, y1, x2, y2);

% Pick intersection closest to y = x
[~, mid] = min(abs(xi - yi));
P2x = xi(mid); 
P2y = yi(mid);

fprintf('Metastable: (%.4f, %.4f)\n', P2x, P2y);

%% Interpolators
[x1s, ix1] = sort(x1); y1s = y1(ix1);
[x2s, ix2] = sort(x2); y2s = y2(ix2);

c1_y = @(x) interp1(x1s, y1s, x, 'linear', NaN); % upper curve
c2_y = @(x) interp1(x2s, y2s, x, 'linear', NaN); % lower curve

%% =========================
%% LOBE 1 (LEFT SIDE)
%% =========================
best_s1 = 0;
best_sq1 = [0 0 0 0];

%idx1 = find(x1 < P2x); % curve1 (top)
%idx2 = find(x2 < P2x); % curve2 (bottom)

tol = 1e-3;
for bx = linspace(0, P2x, 300)
    for s = linspace(0.001, 0.5, 300)

        by = min(c1_y(bx), c2_y(bx)); % bottom inside lobe
        if isnan(by), continue; end

        top = by + s;

        % Check BOTH top corners
        yL = max(c1_y(bx), c2_y(bx));
        yR = max(c1_y(bx+s), c2_y(bx+s));

        if top > yL || top > yR
            continue;
        end

        % Check bottom-right stays inside
        y_lower_R = min(c1_y(bx+s), c2_y(bx+s));
        if by < y_lower_R
            continue;
        end

        if s > best_s1
            best_s1 = s;
            best_sq1 = [bx, by, s, s];
        end
    end
end

%{
for i = idx1'
    x_top = x1(i);
    y_top = y1(i);

    for j = idx2'
        x_bot = x2(j);
        y_bot = y2(j);

        s = min(y_top - y_bot, x_bot - x_top);
        if s <= 0, continue; end
        
        % Square corners
        xs = [x_top, x_top+s, x_top, x_top+s];
        ys = [y_bot, y_bot, y_bot+s, y_bot+s];

        valid = true;
        
        for k = 1:4
            xk = xs(k); 
            yk = ys(k);
        
            y1k = c1_y(xk);
            y2k = c2_y(xk);
        
            if isnan(y1k) || isnan(y2k)
                valid = false; break;
            end
        
            y_upper = max(y1k, y2k);
            y_lower = min(y1k, y2k);
        
            if yk > (y_upper + tol) || yk < (y_lower - tol)
                valid = false; break;
            end
        end

        if valid && s > best_s1
            best_s1 = s;
            best_sq1 = [x_top, y_bot, s, s];
        end
    end
end
%}

%% =========================
%% LOBE 2 (RIGHT SIDE)
%% =========================
best_s2 = 0;
best_sq2 = [0 0 0 0];

idx1 = find(x1 > P2x); % curve1
idx2 = find(x2 > P2x); % curve2

for i = idx1'
    x_bot = x1(i);
    y_bot = y1(i);

    for j = idx2'
        x_top = x2(j);
        y_top = y2(j);

        s = min(y_top - y_bot, x_top - x_bot);
        if s <= 0, continue; end

        % Square corners
        xs = [x_bot, x_bot+s, x_bot, x_bot+s];
        ys = [y_bot, y_bot, y_bot+s, y_bot+s];

        tol = 1e-3;
        valid = true;
        
        for k = 1:4
            xk = xs(k); 
            yk = ys(k);
        
            y1k = c1_y(xk);
            y2k = c2_y(xk);
        
            if isnan(y1k) || isnan(y2k)
                valid = false; break;
            end
        
            y_upper = max(y1k, y2k);
            y_lower = min(y1k, y2k);
        
            if yk > y_upper + tol || yk < y_lower - tol
                valid = false; break;
            end
        end

        if valid && s > best_s2
            best_s2 = s;
            best_sq2 = [x_bot, y_bot, s, s];
        end
    end
end

%% SNM
SNM = min(best_s1, best_s2);

fprintf('Lobe1 SNM = %.4f V\n', best_s1);
fprintf('Lobe2 SNM = %.4f V\n', best_s2);
fprintf('Final SNM = %.4f V\n', SNM);

%% =========================
%% Plot
%% =========================
figure; set(gcf,'Color','k');
ax = axes; 
ax.Color='k'; ax.XColor='w'; ax.YColor='w'; ax.GridColor='w';
hold on; grid on; axis equal;

plot(x1,y1,'g','LineWidth',2);
plot(x2,y2,'r','LineWidth',2);

lim = max([x1; x2])*1.05;
plot([0 lim],[0 lim],'w--','LineWidth',1.5);

% Squares
rectangle('Position',best_sq1,'EdgeColor','y','LineWidth',2.5);
rectangle('Position',best_sq2,'EdgeColor','c','LineWidth',2.5);

% Intersections
plot(xi, yi, 'ro', 'MarkerSize',8,'LineWidth',2);

xlabel('Left node','Color','w');
ylabel('Right node','Color','w');
title(sprintf('SNM = %.4f V', SNM),'Color','w');

legend('Vq vs Vqb','Vqb vs Vq','y=x','Lobe1 sq','Lobe2 sq',...
       'TextColor','w','Color','k');
