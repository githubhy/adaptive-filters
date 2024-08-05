function [y, e] = dfe_rls(x, d, M, N, lambda, delta)
    % DFE-RLS ISI cancellation adaptation algorithm for QPSK
    %
    % Inputs:
    % x: input signal (QPSK symbols)
    % d: desired signal (QPSK symbols)
    % M: number of feedforward taps
    % N: number of feedback taps
    % lambda: forgetting factor (typically close to 1, e.g., 0.99)
    % delta: initial value of P (small positive constant)
    %
    % Outputs:
    % y: equalized output signal
    % e: error signal
    %
    % Ref: TABLE 10.1 Summary of the RLS Algorithm in Adaptive Filter
    % Theory

    L = length(x);
    
    % Initialize filter coefficients
    w = zeros(M+N, 1);
    
    % Initialize inverse correlation matrix
    P = (1/delta) * eye(M+N);
    
    % Initialize output and error signals
    y = zeros(L, 1);
    e = zeros(L, 1);
    
    % Buffer for previous decisions
    decision_buffer = zeros(N, 1);
    
    for n = M:L
        % Construct input vector
        x_n = [x(n:-1:n-M+1); decision_buffer];
        
        % Calculate output
        y(n) = w' * x_n;
        
        % Make decision for QPSK
        decision = qpsk_decision(y(n));
        
        % Calculate error
        e(n) = d(n) - y(n);
        
        % Compute gain vector
        k = (P * x_n) / (lambda + x_n' * P * x_n);
        
        % Update inverse correlation matrix
        P = (P - k * x_n' * P) / lambda;
        
        % Update filter coefficients
        w = w + k * conj(e(n));
        
        % Update decision buffer
        decision_buffer = [decision; decision_buffer(1:end-1)];
    end
end