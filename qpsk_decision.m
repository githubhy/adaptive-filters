function decision = qpsk_decision(symbol)
    % QPSK decision function
    if real(symbol) >= 0 && imag(symbol) >= 0
        decision = (1 + 1i) / sqrt(2);
    elseif real(symbol) < 0 && imag(symbol) >= 0
        decision = (-1 + 1i) / sqrt(2);
    elseif real(symbol) < 0 && imag(symbol) < 0
        decision = (-1 - 1i) / sqrt(2);
    else
        decision = (1 - 1i) / sqrt(2);
    end
end