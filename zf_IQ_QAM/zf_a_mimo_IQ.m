clear all;
close all;
clc
% 2发2收的MIMO，发射采用 Alamouti STBC ，zf预编码后在模拟域调制输出
% 不包含信道估计
jay = sqrt(-1);
fc = 10e3;
fs = 50e3;
QAM_order = 16;
k = log2(QAM_order);
symble_per_time = 1000;
sample_per_symbol = 100;    % fs = smple_rate * sample_per_symbol
deltaT = 1 / fs;

symbol_rate = fs / sample_per_symbol;
SigTim = symble_per_time / symbol_rate;
DataTim = 0 : deltaT : SigTim - deltaT;
CarrierSig = exp(jay * 2 * pi * fc * DataTim);
deCos = cos(2 * pi * fc * (0 : sample_per_symbol - 1) / fs);
deSin =  -sin(2 * pi * fc * (0 : sample_per_symbol - 1) / fs);

randi_dec = randi(QAM_order, symble_per_time, 1) - 1;
randi_bit = de2bi(randi_dec, k);    % 误码率统计用
info = qammod(randi_dec, QAM_order);

%% - - - Channel- - - %%
H_base = [1, 0; 0, 1];
H_fliplr = conj(fliplr(H_base));
H_fliplr(:, 2) = -H_fliplr(:, 2);
H = [H_base; H_fliplr];
inv_H = (H' * H) \ H';

%% - - - STBC - - - %%
stbc_temp = kron(reshape(info, 2, []), ones(1, 2));
stbc_temp(1 , 2 : 2 : end) = conj(stbc_temp(1 , 2 : 2 : end));
stbc_temp(2 , 2 : 2 : end) = -conj(stbc_temp(2 , 2 : 2 : end));
stbc_temp(:, 2: 2: end) = flip(stbc_temp(:, 2: 2: end));

%% - - - ZF - - - %%
coded_info = stbc_temp;

%% - - - IQ Mod - - - %%
CarrierSig_reshape = reshape(CarrierSig, sample_per_symbol, []);
for counter = 1 : length(CarrierSig_reshape)
   IQ_zf_mod_temp(1, :, counter) = coded_info(1, counter) * CarrierSig_reshape(:, counter);
   IQ_zf_mod_temp(2, :, counter) = coded_info(2, counter) * CarrierSig_reshape(:, counter);
end

IQ_zf_mod(1, :) = real(reshape(IQ_zf_mod_temp(1, :, :), 1, []));
IQ_zf_mod(2, :) = real(reshape(IQ_zf_mod_temp(2, :, :), 1, []));

%% - - - RX - - - %%
for counter_i = 1 : length(IQ_zf_mod)
    rx_temp(:, counter_i) = H_base * IQ_zf_mod(:, counter_i);
end

%% - - - IQ Demod - - - %%
rx_data.one = reshape(rx_temp(1, :), 100, []);
rx_data.two = reshape(rx_temp(2, :), 100, []);
rx_zf_info(1, :) = demodIQ(deCos, deSin, rx_data.one, sample_per_symbol);
rx_zf_info(2, :) = demodIQ(deCos, deSin, rx_data.two, sample_per_symbol);

%% - - - STBC Estimate - - - %%
% 通过信道补偿估计STBC的值
% for counter_i = 1 : length(rx_zf_info)
% rx_est_stbc(:, counter_i) = inv_H * rx_zf_info(:, counter_i);
% end
rx_dt = rx_zf_info(:, 1 : 2 : end);
rx_dt = [rx_dt; conj(rx_zf_info(:, 2 : 2 : end))];

stbc_est = inv_H * rx_dt;
stbc_est(2, :) = stbc_est(2, :);

%% - - - QAM Demod - - - %%
rx_zf_qam_demod(1, :) = qamdemod(stbc_est(1, :), QAM_order);
rx_zf_qam_demod(2, :) = qamdemod(stbc_est(2, :), QAM_order);
rx_decode = zeros(1, 2 * length(rx_zf_qam_demod));
rx_decode(1 : 2 :end) = rx_zf_qam_demod(1, :);
rx_decode(2 : 2 :end) = rx_zf_qam_demod(2, :);
rx_in_bit = de2bi(rx_decode);

%% - - - Error - - - %%
[numErrors, ber] = biterr(rx_in_bit(:), randi_bit(:));

%% - - - Plot - - - %%
% figure(1)
% plot(real(stbc_temp), imag(stbc_temp), 'co');
% hold on;
% plot(real(rx_est_stbc), imag(rx_est_stbc), 'r*');
% hold on;
% plot(real(rx_zf_info), imag(rx_zf_info), 'd')
