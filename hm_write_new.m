clc; clear; close all;

%% Load data
data1 = readtable('read2.txt',  'VariableNamingRule','preserve');
data2 = readtable('write2.txt', 'VariableNamingRule','preserve');

Vq1=data1{:,2}; Vqb1=data1{:,3};
Vq2=data2{:,2}; Vqb2=data2{:,3};

x1=Vq1;  y1=Vqb1;   % green curve (read)
x2=Vqb2; y2=Vq2;    % red curve   (write)

%% Build x(y) interpolators
[y1s,iy1]=sort(y1); x1s=x1(iy1); [y1s,u1]=unique(y1s); x1s=x1s(u1);
[y2s,iy2]=sort(y2); x2s=x2(iy2); [y2s,u2]=unique(y2s); x2s=x2s(u2);

c_green = @(y) interp1(y1s, x1s, y, 'linear', NaN);  % right boundary
c_red   = @(y) interp1(y2s, x2s, y, 'linear', NaN);  % left boundary

%% Dense y grid
y_lo = max(y1s(1),   y2s(1));
y_hi = min(y1s(end), y2s(end));
N    = 5000;
yg   = linspace(y_lo, y_hi, N)';
xG   = c_green(yg);
xR   = c_red(yg);

%% Find largest square where:
%%   available_width = min(xG) - max(xR)  over [yb, yb+s]  >= s
best_s  = 0;
best_sq = [0 0 0 0];

for i = 1:N
    yb   = yg(i);
    s_lo = 0;
    s_hi = y_hi - yb;

    for iter = 1:80
        s_mid = (s_lo + s_hi) / 2;
        mask  = yg >= yb & yg <= yb + s_mid;
        if ~any(mask), s_hi = s_mid; continue; end

        % KEY FIX: available width = min(green) - max(red) over full height range
        avail = min(xG(mask)) - max(xR(mask));

        if avail >= s_mid
            s_lo = s_mid;   % fits — try bigger
        else
            s_hi = s_mid;   % doesn't fit — shrink
        end
    end

    s = s_lo;
    if s > best_s
        mask    = yg >= yb & yg <= yb + s;
        xl      = max(xR(mask));    % push square right of tightest red point
        best_s  = s;
        best_sq = [xl, yb, s, s];
    end
end

SNM = best_s;
fprintf('SNM = %.4f V\n', SNM);

%% Plot
figure; set(gcf,'Color','k');
ax=axes; ax.Color='k'; ax.XColor='w'; ax.YColor='w'; ax.GridColor='w';
hold on; grid on; axis equal;

plot(x1,y1,'g','LineWidth',2);
plot(x2,y2,'r','LineWidth',2);
rectangle('Position',best_sq,'EdgeColor','y','LineWidth',2.5);

xlabel('Left node','Color','w');
ylabel('Right node','Color','w');
title(sprintf('SNM = %.4f V', SNM),'Color','w');
legend('Read (green)','Write (red)','SNM square','TextColor','w','Color','k');