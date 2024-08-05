function decision = qpsk_decision_fixed(symbol)
    % Fixed-point QPSK decision function
    ft_decision = fixdt(1, 16, 14);  % Q1.14 for decision
    
    real_part = real(symbol);
    imag_part = imag(symbol);
    
    if real_part >= 0 && imag_part >= 0
        decision = fi((1 + 1i) / sqrt(2), ft_decision);
    elseif real_part < 0 && imag_part >= 0
        decision = fi((-1 + 1i) / sqrt(2), ft_decision);
    elseif real_part < 0 && imag_part < 0
        decision = fi((-1 - 1i) / sqrt(2), ft_decision);
    else
        decision = fi((1 - 1i) / sqrt(2), ft_decision);
    end
end