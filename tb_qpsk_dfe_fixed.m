clearvars;
close all;

% Parameters
num_symbols = 100;  % Number of symbols
mu_ff = 0.01;  % Step size for LMS adaptation
mu_fb = 0.01;
lambda = 0.99;  % Forgetting factor for RLS
delta = 0.1;  % Initial value of P for RLS

% SNR range
snr_range = -10:2:20;
num_snr_points = length(snr_range);

% M and N values
M = 16;
N = 16;

% Preallocate arrays for BER results
ber_no_eq = zeros(1, num_snr_points);
ber_lms_float = zeros(1, num_snr_points);
ber_rls_float = zeros(1, num_snr_points);
ber_lms_fixed = zeros(1, num_snr_points);
ber_rls_fixed = zeros(1, num_snr_points);

% Generate random TX bits with QPSK modulation
tx_bits = randi([0 1], 2*num_symbols, 1);
tx_symbols = qpsk_modulate(tx_bits);

% Define and normalize multipath channel
channel = [1 0.5 0.3 0 0.2 0 0 0.1];
channel = channel / norm(channel);

% Apply multipath channel
rx_symbols = conv(tx_symbols, channel, 'same');

% Parallel processing for SNR range
for i = 1:num_snr_points
    snr_db = snr_range(i);
    
    % Add noise
    rx_symbols_noisy = awgn(rx_symbols, snr_db, 'measured');

    % Calculate BER without equalization
    rx_bits_no_eq = qpsk_demodulate(rx_symbols_noisy);
    ber_no_eq(i) = sum(rx_bits_no_eq ~= tx_bits) / length(tx_bits);

    % Calculate BER with floating-point DFE-LMS equalization
    [eq_symbols_lms_float, ~] = dfe_lms(rx_symbols_noisy, tx_symbols, M, N, mu_ff, mu_fb);
    rx_bits_lms_float = qpsk_demodulate(eq_symbols_lms_float);
    ber_lms_float(i) = sum(rx_bits_lms_float ~= tx_bits) / length(tx_bits);

    % Calculate BER with floating-point DFE-RLS equalization
    [eq_symbols_rls_float, ~] = dfe_rls(rx_symbols_noisy, tx_symbols, M, N, lambda, delta);
    rx_bits_rls_float = qpsk_demodulate(eq_symbols_rls_float);
    ber_rls_float(i) = sum(rx_bits_rls_float ~= tx_bits) / length(tx_bits);

    % Calculate BER with fixed-point DFE-LMS equalization
    [eq_symbols_lms_fixed, ~] = dfe_lms_fixed(rx_symbols_noisy, tx_symbols, M, N, mu_ff, mu_fb);
    rx_bits_lms_fixed = qpsk_demodulate(double(eq_symbols_lms_fixed));
    ber_lms_fixed(i) = sum(rx_bits_lms_fixed ~= tx_bits) / length(tx_bits);

    % Calculate BER with fixed-point DFE-RLS equalization
    [eq_symbols_rls_fixed, ~] = dfe_rls_fixed(rx_symbols_noisy, tx_symbols, M, N, lambda, delta);
    rx_bits_rls_fixed = qpsk_demodulate(double(eq_symbols_rls_fixed));
    ber_rls_fixed(i) = sum(rx_bits_rls_fixed ~= tx_bits) / length(tx_bits);

    fprintf('Completed simulation for SNR(dB)=%f.\n', snr_db);
end

% Plot BER vs SNR for all configurations
figure;
semilogy(snr_range, ber_no_eq, 'k-', 'LineWidth', 2);
hold on;
semilogy(snr_range, ber_lms_float, 'ro-', 'LineWidth', 1.5);
semilogy(snr_range, ber_rls_float, 'bo-', 'LineWidth', 1.5);
semilogy(snr_range, ber_lms_fixed, 'rs--', 'LineWidth', 1.5);
semilogy(snr_range, ber_rls_fixed, 'bs--', 'LineWidth', 1.5);

grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for QPSK with Floating-point and Fixed-point DFE Equalization');
legend('No Equalization', 'LMS Float', 'RLS Float', 'LMS Fixed', 'RLS Fixed', 'Location', 'southwest');
ylim([1e-5, 1]);

% Print best performance for each method
[min_ber_lms_float, best_snr_lms_float] = min(ber_lms_float);
[min_ber_rls_float, best_snr_rls_float] = min(ber_rls_float);
[min_ber_lms_fixed, best_snr_lms_fixed] = min(ber_lms_fixed);
[min_ber_rls_fixed, best_snr_rls_fixed] = min(ber_rls_fixed);

fprintf('Best LMS Float: SNR=%d dB, BER=%e\n', snr_range(best_snr_lms_float), min_ber_lms_float);
fprintf('Best RLS Float: SNR=%d dB, BER=%e\n', snr_range(best_snr_rls_float), min_ber_rls_float);
fprintf('Best LMS Fixed: SNR=%d dB, BER=%e\n', snr_range(best_snr_lms_fixed), min_ber_lms_fixed);
fprintf('Best RLS Fixed: SNR=%d dB, BER=%e\n', snr_range(best_snr_rls_fixed), min_ber_rls_fixed);
