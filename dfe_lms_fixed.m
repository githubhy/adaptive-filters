function [y, e] = dfe_lms_fixed(x, d, M, N, mu_ff, mu_fb)
    % Fixed-point DFE-LMS ISI cancellation adaptation algorithm for QPSK
    
    % Define fixed-point data types
    ft_input = fixdt(1, 16, 14);  % Q1.14 for input
    ft_coeff = fixdt(1, 16, 14);  % Q1.14 for coefficients
    ft_acc = fixdt(1, 32, 30);    % Q1.30 for accumulator
    
    L = length(x);
    
    % Initialize filter coefficients
    w_ff = fi(zeros(M, 1), ft_coeff);
    w_fb = fi(zeros(N, 1), ft_coeff);
    
    % Initialize output and error signals
    y = fi(zeros(L, 1), ft_input);
    e = fi(zeros(L, 1), ft_input);
    
    % Buffer for previous decisions
    decision_buffer = fi(zeros(N, 1), ft_input);
    
    for n = M:L
        % Feedforward section
        x_ff = fi(x(n:-1:n-M+1), ft_input);
        
        % Feedback section
        x_fb = decision_buffer;
        
        % Calculate output
        y(n) = fi(sum(w_ff .* x_ff) + sum(w_fb .* x_fb), ft_acc);
        
        % Make decision for QPSK
        decision = qpsk_decision_fixed(y(n));
        
        % Calculate error
        e(n) = fi(d(n) - y(n), ft_input);
        
        % Update filter coefficients
        w_ff = fi(w_ff + fi(mu_ff * e(n) * conj(x_ff), ft_coeff), ft_coeff);
        w_fb = fi(w_fb + fi(mu_fb * e(n) * conj(x_fb), ft_coeff), ft_coeff);
        
        % Update decision buffer
        decision_buffer = [decision; decision_buffer(1:end-1)];
    end
end