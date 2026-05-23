% nn_mac_noise.m
% -----------------------------------------------------------------------------
% MNIST 2-layer NN pipeline (from document.xml) with Gaussian noise injected
% into the MAC stage to model the analog mismatch characterized by
% sig_delta_mac_plot_final.m. That companion script produces
% `sigma_over_delta`, the RMS MAC error expressed in LSBs of the 5-bit spice
% ADC. Here we sweep that value, run NUM_TRIALS Monte-Carlo trials per sigma,
% and report mean +/- std test accuracy on MNIST.
%
% Noise model: for each chunk a = 1..8, MAC{a} is a 100x10000 matrix of
% pre-quantization MAC results. We add independent N(0, sigma_NN^2) noise to
% every entry BEFORE the clip + fi(_,1,5,0) quantization. sigma_NN is the
% noise std expressed in NN-LSB units; since the NN's 5-bit signed quantizer
% has LSB = 1 MAC count, sigma_NN = sigma_over_delta * LSB_scale.
% Set LSB_scale = 1 for NN-LSB scaling, or 8 if you instead want to interpret
% sigma_over_delta in spice-ADC LSBs (which are ~8 MAC counts wide given the
% alpha = 13.5 clip and 5-bit signed range).
% -----------------------------------------------------------------------------

clear;  clc;

% -------------------- Tunable parameters --------------------
sigma_over_delta_list = [0, 1/sqrt(12), 0.4, 0.6, 0.8, 1.0, 1.5, 2.0];
NUM_TRIALS            = 5;
LSB_scale             = 1;     % MAC counts per NN LSB; set to 8 for spice-LSB scaling
RNG_SEED              = 20260523;

rng(RNG_SEED);

% -------------------- Load data + weights --------------------
load('weightsAndBiases.mat', 'W1', 'b1', 'W2', 'b2');
images = loadMNISTImages('t10k-images');
labels = loadMNISTLabels('t10k-labels')';

% -------------------- Quantize, pad, split (once) --------------------
w = 6;  wf = 2;  x = 2;  xf = 2;          % bitwidths
W1_Q = fi(W1, 1, w, wf);
X1_Q = fi(images, 0, x, xf);

% Pad W1 to 1024 cols, X to 1024 rows; split into 8 chunks of 128
W_Q     = [W1_Q fi(zeros(size(W1_Q,1), 1024 - size(W1_Q,2)), 1, w, wf)];
W_split = mat2cell(W_Q, size(W_Q,1), repmat(128,1,8));
X_Q     = [X1_Q; fi(zeros(1024 - size(X1_Q,1), size(X1_Q,2)), 0, x, xf)];
X_split = mat2cell(X_Q, repmat(128,1,8), size(X_Q,2));

% -------------------- Compute noise-free MAC{a} once --------------------
MAC = cell(1, 8);
for a = 1:8
    MAC{a} = W_split{a}.data * X_split{a}.data;   % 100 x 10000
end

% Clip + quantize parameters
q     = 5;
alpha = 13.5;
[nRows, nCols] = size(MAC{1});                    % 100 x 10000

% -------------------- Sweep sigma, Monte-Carlo trials --------------------
nSig       = numel(sigma_over_delta_list);
mean_acc   = zeros(1, nSig);
std_acc    = zeros(1, nSig);
all_acc    = zeros(nSig, NUM_TRIALS);

fprintf('Running %d sigma values x %d trials each...\n', nSig, NUM_TRIALS);

for s = 1:nSig
    sigma_over_delta = sigma_over_delta_list(s);
    sigma_NN         = sigma_over_delta * LSB_scale;   % std-dev in MAC counts

    for t = 1:NUM_TRIALS
        finalMAC = zeros(nRows, nCols);
        for a = 1:8
            if sigma_NN > 0
                MAC_noisy = MAC{a} + randn(nRows, nCols) * sigma_NN;
            else
                MAC_noisy = MAC{a};                 % noise-free baseline
            end
            MAC_clipped_a = min(max(MAC_noisy, -alpha), alpha);
            MAC_Q_a       = fi(MAC_clipped_a, 1, q, 0);
            finalMAC      = finalMAC + MAC_Q_a.data;
        end
        finalMAC = finalMAC + b1;
        a1       = max(0, finalMAC);

        % Layer 2
        z2 = W2 * a1 + b2;
        [~, guesses] = max(z2, [], 1);
        guesses      = guesses - 1;
        Acc          = sum(labels == guesses) / numel(labels) * 100;

        all_acc(s, t) = Acc;
        fprintf('  sigma/delta = %6.4f  trial %d/%d :  Acc = %6.3f%%\n', ...
            sigma_over_delta, t, NUM_TRIALS, Acc);
    end

    mean_acc(s) = mean(all_acc(s, :));
    std_acc(s)  = std (all_acc(s, :));
end

% -------------------- Summary table --------------------
fprintf('\n=================== Summary ===================\n');
fprintf(' sigma/delta |  mean Acc (%%) |  std Acc (%%)\n');
fprintf(' ------------+---------------+--------------\n');
for s = 1:nSig
    fprintf('  %9.4f  |   %9.4f   |  %9.4f\n', ...
        sigma_over_delta_list(s), mean_acc(s), std_acc(s));
end
fprintf('===============================================\n');

% Noise-free baseline (sigma_over_delta == 0)
baseline_idx = find(sigma_over_delta_list == 0, 1, 'first');
if isempty(baseline_idx)
    baseline_acc = NaN;
else
    baseline_acc = mean_acc(baseline_idx);
end

% -------------------- Plot --------------------
figure;
errorbar(sigma_over_delta_list, mean_acc, std_acc, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'auto');
hold on;
if ~isnan(baseline_acc)
    yline(baseline_acc, '--', sprintf('noise-free = %.2f%%', baseline_acc), ...
        'LineWidth', 1.2);
end
hold off;
grid on;
xlabel('\sigma / \Delta  (MAC-error RMS in NN LSBs)');
ylabel('MNIST test accuracy (%)');
title(sprintf('NN accuracy vs MAC noise  (%d trials/point, LSB scale = %g)', ...
    NUM_TRIALS, LSB_scale));
