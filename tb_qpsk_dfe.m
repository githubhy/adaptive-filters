clearvars;
close all;

% Parameters
num_symbols = 100000;  % Number of symbols
mu_ff = 0.01;  % Step size for LMS adaptation
mu_fb = 0.01;
lambda = 0.99;  % Forgetting factor for RLS
delta = 0.1;  % Initial value of P for RLS

% SNR range
snr_range = -10:2:20;
num_snr_points = length(snr_range);

% M and N ranges (simplified for comparison)
M_range = [8, 16];
N_range = [8, 16];

% Preallocate arrays for BER results
ber_no_eq = zeros(1, num_snr_points);
ber_lms = zeros(length(M_range), length(N_range), num_snr_points);
ber_rls = zeros(length(M_range), length(N_range), num_snr_points);

% Generate random TX bits with QPSK modulation
tx_bits = randi([0 1], 2*num_symbols, 1);
tx_symbols = qpsk_modulate(tx_bits);

% Define and normalize multipath channel
channel = [1 0.5 0.3 0 0.2 0 0 0.1];
channel = channel / norm(channel);

% Apply multipath channel
rx_symbols = conv(tx_symbols, channel, 'same');

% Main loop for different M and N combinations
for m = 1:length(M_range)
    for n = 1:length(N_range)
        M = M_range(m);
        N = N_range(n);
        
        % Parallel processing for SNR range
        parfor i = 1:num_snr_points
            snr_db = snr_range(i);
            
            % Add noise
            rx_symbols_noisy = awgn(rx_symbols, snr_db, 'measured');

            % Calculate BER without equalization (only for the first M,N combination)
            if m == 1 && n == 1
                rx_bits_no_eq = qpsk_demodulate(rx_symbols_noisy);
                ber_no_eq(i) = sum(rx_bits_no_eq ~= tx_bits) / length(tx_bits);
            end

            % Calculate BER with DFE-LMS equalization
            [eq_symbols_lms, ~] = dfe_lms(rx_symbols_noisy, tx_symbols, M, N, mu_ff, mu_fb);
            rx_bits_lms = qpsk_demodulate(eq_symbols_lms);
            ber_lms(m, n, i) = sum(rx_bits_lms ~= tx_bits) / length(tx_bits);

            % Calculate BER with DFE-RLS equalization
            [eq_symbols_rls, ~] = dfe_rls(rx_symbols_noisy, tx_symbols, M, N, lambda, delta);
            rx_bits_rls = qpsk_demodulate(eq_symbols_rls);
            ber_rls(m, n, i) = sum(rx_bits_rls ~= tx_bits) / length(tx_bits);
        end
        
        fprintf('Completed simulation for M=%d, N=%d\n', M, N);
    end
end

% Plot BER vs SNR for all configurations
figure;
semilogy(snr_range, ber_no_eq, 'k-', 'LineWidth', 2);
hold on;

colors = ['r', 'g', 'b', 'c'];
markers = ['o', 's', 'd', '^'];

legend_str = {'Without Equalization'};

for m = 1:length(M_range)
    for n = 1:length(N_range)
        color_idx = (m-1)*length(N_range) + n;
        
        semilogy(snr_range, squeeze(ber_lms(m, n, :)), ...
            [colors(color_idx), markers(1), '-'], 'LineWidth', 1.5);
        semilogy(snr_range, squeeze(ber_rls(m, n, :)), ...
            [colors(color_idx), markers(2), '-'], 'LineWidth', 1.5);
        
        legend_str{end+1} = sprintf('LMS M=%d, N=%d', M_range(m), N_range(n));
        legend_str{end+1} = sprintf('RLS M=%d, N=%d', M_range(m), N_range(n));
    end
end

grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for QPSK with DFE-LMS and DFE-RLS Equalization');
legend(legend_str, 'Location', 'southwest');
ylim([1e-5, 1]);

% Find best configurations for LMS and RLS
[min_ber_lms, idx_lms] = min(ber_lms(:));
[best_m_lms, best_n_lms, best_snr_lms] = ind2sub(size(ber_lms), idx_lms);
[min_ber_rls, idx_rls] = min(ber_rls(:));
[best_m_rls, best_n_rls, best_snr_rls] = ind2sub(size(ber_rls), idx_rls);

fprintf('Best LMS configuration: M=%d, N=%d at SNR=%d dB (BER=%e)\n', ...
    M_range(best_m_lms), N_range(best_n_lms), snr_range(best_snr_lms), min_ber_lms);
fprintf('Best RLS configuration: M=%d, N=%d at SNR=%d dB (BER=%e)\n', ...
    M_range(best_m_rls), N_range(best_n_rls), snr_range(best_snr_rls), min_ber_rls);

%{
% Plot constellation diagrams
figure;
subplot(2,1,1);
scatter(real(rx_symbols_noisy), imag(rx_symbols_noisy), 'r.'); hold on;
scatter(real(tx_symbols), imag(tx_symbols), 'bo');
title('Transmitted and Received Constellations (No Equalization)');
xlabel('In-phase'); ylabel('Quadrature');
axis([-2 2 -2 2]);
legend('Transmitted', 'Received (No EQ)');

subplot(2,1,2);
scatter(real(eq_symbols), imag(eq_symbols), 'r.'); hold on;
scatter(real(tx_symbols), imag(tx_symbols), 'bo');
title('Transmitted and Received Constellations (With DFE-LMS Equalization)');
xlabel('In-phase'); ylabel('Quadrature');
axis([-2 2 -2 2]);
legend('Transmitted', 'Received (With EQ)');

% Plot LMS convergence
figure;
mse = abs(error).^2;
semilogy(movmean(mse, 100));  % Apply moving average for smoother curve
title('LMS Convergence');
xlabel('Symbol Index');
ylabel('Mean Squared Error');
grid on;
%}
