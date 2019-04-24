function [y, len] = channel(x, snr_in_db)
SNR = exp(snr_in_db * log(10) / 10);
sigma = 1 / sqrt(2 * SNR);

actual_isi = [0.05 -0.063 0.088 -0.126 -0.25 0.9047 0.25 0 0.126 0.038 0.088];
len_actual_isi = (length(actual_isi) - 1) / 2;
len = len_actual_isi;
y = conv(actual_isi, x);

for i = 1 : 2 : size(y, 2)
    [noise(i) noise(i+1)] = gngauss(sigma);
end

y = y + noise;
end