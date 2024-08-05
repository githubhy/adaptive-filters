function [y, e] = dfe_rls_fixed(x, d, M, N, lambda, delta)
    % Fixed-point DFE-RLS ISI cancellation adaptation algorithm for QPSK
    
    % Define fixed-point data types
    ft_input = fixdt(1, 16, 14);  % Q1.14 for input
    ft_coeff = fixdt(1, 16, 14);  % Q1.14 for coefficients
    ft_acc = fixdt(1, 32, 30);    % Q1.30 for accumulator
    
    L = length(x);
    
    % Initialize filter coefficients
    w = fi(zeros(M+N, 1), ft_coeff);
    
    % Initialize inverse correlation matrix
    P = fi((1/delta) * eye(M+N), ft_acc);
    
    % Initialize output and error signals
    y = fi(zeros(L, 1), ft_input);
    e = fi(zeros(L, 1), ft_input);
    
    % Buffer for previous decisions
    decision_buffer = fi(zeros(N, 1), ft_input);
    
    for n = M:L
        % Construct input vector
        x_n = fi([x(n:-1:n-M+1); decision_buffer], ft_input);
        
        % Calculate output
        y(n) = fi(w' * x_n, ft_acc);
        
        % Make decision for QPSK
        decision = qpsk_decision_fixed(y(n));
        
        % Calculate error
        e(n) = fi(d(n) - y(n), ft_input);
        
        % Compute gain vector
        dno = lambda + x_n' * P * x_n;
        den = fi(fi(1, ft_acc) / dno, ft_acc);
        k = fi(P * x_n * den, ft_acc);
        
        % Update inverse correlation matrix
        P = fi((P - fi(k * (x_n' * P), ft_acc)) / lambda, ft_acc);
        
        % Update filter coefficients
        w = fi(w + fi(k * conj(e(n)), ft_coeff), ft_coeff);
        
        % Update decision buffer
        decision_buffer = [decision; decision_buffer(1:end-1)];
    end
end
