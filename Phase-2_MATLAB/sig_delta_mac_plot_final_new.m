clear; clc;

% -------------------------------------------------------------------------
% MAC Quantization Analysis
%
% Set STEP_SELECT below to choose which step to read from mac_neg / mac_pos
%   1 = step 1
%   2 = step 2
% -------------------------------------------------------------------------

STEP_SELECT = 2;   % <-- change this to 1 or 2

% -------------------------------------------------------------------------

bits       = 5;
num_levels = 2^bits;
num_mac    = 128;

%% -----------------------------------------------------------------------
%  PARSER 1: old format  "v_sampleN: v(v_out)=<val> at <time>"
% -----------------------------------------------------------------------
function v = read_simple(filename, num_mac)
    fid = fopen(filename, 'r');
    if fid == -1, error('Cannot open: %s', filename); end
    v = zeros(1, num_mac); idx = 0;
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if ischar(line)
            tok = regexp(line, 'v\(v_out\)=([\d.\-e+]+)', 'tokens');
            if ~isempty(tok)
                idx = idx + 1;
                v(idx) = str2double(tok{1}{1});
            end
        end
    end
    fclose(fid);
    v = v(1:idx);
end

%% -----------------------------------------------------------------------
%  PARSER 2: new format — selectable step
%  "Measurement: V_sampleN" ... "1   <volt>   <time>"
%                              "2   <volt>   <time>"
% -----------------------------------------------------------------------
function v = read_step(filename, num_mac, step_sel)
    fid = fopen(filename, 'r');
    if fid == -1, error('Cannot open: %s', filename); end
    v = zeros(1, num_mac); idx = 0;
    in_meas = false;
    step_str = num2str(step_sel);
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if ischar(line)
            if startsWith(line, 'Measurement:')
                in_meas = true;
            elseif in_meas && startsWith(line, step_str)
                parts = sscanf(line, '%d %f %f');
                if numel(parts) >= 2 && parts(1) == step_sel
                    idx = idx + 1;
                    v(idx) = parts(2);
                    in_meas = false;
                end
            end
        end
    end
    fclose(fid);
    v = v(1:idx);
end

%% -----------------------------------------------------------------------
%  1. IDEAL CURVE — from mac_volt / mac_volt_2
% -----------------------------------------------------------------------
v_pos_ideal = read_simple('mac_volt.txt',   num_mac);
v_neg_ideal = read_simple('mac_volt_2.txt', num_mac);
v_neg_ideal = fliplr(v_neg_ideal);

MAC_ideal = [-128:-1,  1:128];
v_ideal   = [v_neg_ideal, v_pos_ideal];

v_min = min(v_ideal);
v_max = max(v_ideal);
delta = (v_max - v_min) / num_levels;

q_idx_ideal   = min(floor((v_ideal - v_min) / delta), num_levels - 1);
v_quant_ideal = v_min + (q_idx_ideal + 0.5) * delta;

%% -----------------------------------------------------------------------
%  2. ACTUAL MEASUREMENTS — mac_neg / mac_pos, selected step
% -----------------------------------------------------------------------
v_neg_act = read_step('mac_neg.txt', num_mac, STEP_SELECT);
v_pos_act = read_step('mac_pos.txt', num_mac, STEP_SELECT);
v_neg_act = fliplr(v_neg_act);

MAC_act = [-128:-1,  1:128];
v_act   = [v_neg_act, v_pos_act];

%% -----------------------------------------------------------------------
%  3. QUANTIZATION ERROR
% -----------------------------------------------------------------------
v_min_act = min(v_act);
v_max_act = max(v_act);
delta_act = (v_max_act - v_min_act) / num_levels;
q_idx_act   = min(floor((v_act - v_min_act) / delta_act), num_levels - 1);
v_quant_act = v_min_act + (q_idx_act + 0.5) * delta_act;

err              = v_act - v_quant_act;
sigma            = sqrt(mean(err.^2));
sigma_over_delta = sigma / delta_act;
ENOB             = bits - log2(sqrt(12) * sigma_over_delta);

fprintf('--- Step %d selected ---\n', STEP_SELECT);
fprintf('  V_min : %.4f V,  V_max : %.4f V\n', v_min, v_max);
fprintf('  Delta (LSB) : %.4f mV\n', delta * 1e3);
fprintf('  Sigma (RMS) : %.4f mV\n', sigma * 1e3);
fprintf('  Sigma/Delta : %.4f  (ideal 1/sqrt(12) = %.4f)\n', ...
        sigma_over_delta, 1/sqrt(12));
fprintf('  ENOB        : %.3f bits  (nominal %d bits)\n', ENOB, bits);

%% -----------------------------------------------------------------------
%  4. PLOT
% -----------------------------------------------------------------------
figure('Name', sprintf('MAC Quantization — Step %d', STEP_SELECT), ...
       'NumberTitle', 'off', 'Color', 'w', 'Position', [100 100 900 500]);

plot(MAC_ideal, v_ideal * 1e3, 'g-', 'LineWidth', 2);
hold on;
plot(MAC_ideal, v_quant_ideal * 1e3, 'r-', 'LineWidth', 1.2);
plot(MAC_act,   v_act * 1e3,   'b--', 'LineWidth', 1.2);
hold off;

xlabel('MAC');
ylabel('V_{out} (mV)');
title(sprintf('Comparison between MAC and Quantized V_{out}  (step %d)', STEP_SELECT));
legend('Ideal V_{out}', 'Quantized (ideal)', ...
       sprintf('Actual V_{out} (step %d)', STEP_SELECT), 'Location', 'best');
grid on;
xlim([-128 128]);

text(-120, v_min*1e3 + 0.08*(v_max - v_min)*1e3, ...
     sprintf('%d-bit \\sigma/\\Delta: %.2f', bits, sigma_over_delta), ...
     'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);
text(-120, v_min*1e3 + 0.02*(v_max - v_min)*1e3, ...
     sprintf('ENOB: %.2f bits', ENOB), ...
     'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);