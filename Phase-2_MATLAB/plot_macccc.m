clear; clc;

% -------------------------------------------------------------------------
% MAC Quantization Analysis
% mac_volt.txt  : MAC =   1 to  128 (128 samples)
% mac_volt_2.txt: MAC =  -1 to -128 (128 samples, stored in reverse order)
% Together they form one continuous transfer curve from -128 to +128
% -------------------------------------------------------------------------

%% --- Configuration ------------------------------------------------------
bits       = 5;
num_levels = 2^bits;      % 32 quantization levels
num_mac    = 128;

%% --- Helper: read one file ----------------------------------------------
function v = read_vout(filename, num_mac)
    fid = fopen(filename, 'r');
    if fid == -1, error('Cannot open file: %s', filename); end
    v = zeros(1, num_mac);
    idx = 0;
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line) && ~isempty(strtrim(line))
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

%% --- Read files ---------------------------------------------------------
v_pos = read_vout('mac_volt.txt',   num_mac);   % MAC 1..128
v_neg = read_vout('mac_volt_2.txt', num_mac);   % MAC -1..-128, flip to -128..-1

% Flip file 2 so it runs from MAC -128 to -1
v_neg = fliplr(v_neg);

% Concatenate: -128..-1 | 1..128  (MAC=0 not present, skip it)
MAC   = [-128:-1,  1:128];
v_out = [v_neg,    v_pos];

%% --- Quantization over full combined range ------------------------------
v_min = min(v_out);
v_max = max(v_out);
delta = (v_max - v_min) / num_levels;

q_idx   = floor((v_out - v_min) / delta);
q_idx   = min(q_idx, num_levels - 1);
v_quant = v_min + (q_idx + 0.5) * delta;

%% --- Sigma / Delta ------------------------------------------------------
err              = v_out - v_quant;
sigma            = sqrt(mean(err.^2));
sigma_over_delta = sigma / delta;

fprintf('--- Quantization Results ---\n');
fprintf('  Voltage range : %.4f V  to  %.4f V\n', v_min, v_max);
fprintf('  Delta (LSB)   : %.4f mV\n', delta * 1e3);
fprintf('  Sigma (RMS)   : %.4f mV\n', sigma * 1e3);
fprintf('  Sigma / Delta : %.4f  (ideal 1/sqrt(12) = %.4f)\n', ...
        sigma_over_delta, 1/sqrt(12));

%% --- Plot ---------------------------------------------------------------
figure('Name', 'MAC Quantization Analysis', 'NumberTitle', 'off', ...
       'Color', 'w', 'Position', [100 100 900 500]);

plot(MAC, v_out   * 1e3, 'g-', 'LineWidth', 2);
hold on;
plot(MAC, v_quant * 1e3, 'r-', 'LineWidth', 1.2);
hold off;

xlabel('MAC');
ylabel('Quantized V_{out}');
title('Comparison between MAC and Quantized V_{out}');
legend('Original V_{out}', 'Quantized V_{out}', 'Location', 'best');
grid on;
xlim([-128 128]);

text(-120, v_min*1e3 + 0.08*(v_max - v_min)*1e3, ...
     sprintf('%d-bit \\sigma/\\Delta: %.2f', bits, sigma_over_delta), ...
     'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);