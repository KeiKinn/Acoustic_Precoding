function demod_data = demodIQ(deCos, deSin, data, sample_per_symbol)
I_info_demod = 2 * (deCos *  data) / sample_per_symbol;
Q_info_demod = 2 * (deSin *  data) / sample_per_symbol;
demod_data = I_info_demod + 1i * Q_info_demod;
end