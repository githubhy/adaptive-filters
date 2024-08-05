function [y, e] = dfe_lms(x, d, M, N, mu_ff, mu_fb)
    % DFE-LMS ISI cancellation adaptation algorithm for QPSK
    %
    % Inputs:
    % x: input signal (QPSK symbols)
    % d: desired signal (QPSK symbols)
    % M: number of feedforward taps
    % N: number of feedback taps
    % mu: step size for LMS adaptation
    %
    % Outputs:
    % y: equalized output signal
    % e: error signal

    L = length(x);
    
    % Initialize filter coefficients
    w_ff = zeros(M, 1);  % Feedforward filter
    w_fb = zeros(N, 1);  % Feedback filter
    
    % Initialize output and error signals
    y = zeros(L, 1);
    e = zeros(L, 1);
    
    % Buffer for previous decisions
    decision_buffer = zeros(N, 1);
    
    for n = M:L
        % Feedforward section
        x_ff = x(n:-1:n-M+1);
        
        % Feedback section
        x_fb = decision_buffer;
        
        % Calculate output
        y(n) = w_ff' * x_ff + w_fb' * x_fb;
        
        % Make decision for QPSK
        decision = qpsk_decision(y(n));
        
        % Calculate error
        e(n) = d(n) - y(n);
        
        % Update filter coefficients
        w_ff = w_ff + mu_ff * x_ff * conj(e(n));
        w_fb = w_fb + mu_fb * x_fb * conj(e(n));
        
        % Update decision buffer
        decision_buffer = [decision; decision_buffer(1:end-1)];
    end
end