function bits = qpsk_demodulate(symbols)
    bits = zeros(2*length(symbols), 1);
    bits(1:2:end) = real(symbols) > 0;
    bits(2:2:end) = imag(symbols) > 0;
end