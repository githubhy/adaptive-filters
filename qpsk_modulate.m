function symbols = qpsk_modulate(bits)
    symbols = (1/sqrt(2)) * ((2*bits(1:2:end)-1) + 1j*(2*bits(2:2:end)-1));
end