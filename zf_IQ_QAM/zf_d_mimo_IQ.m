clear all;
close all;
clc
% 2发2收的MIMO，发射采用 Alamouti STBC ，主要是在数字域中进行理论仿真，暂不添加噪声与预编码
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
H_base = [1, 0.5; 0.8, 0.7];
H_fliplr = conj(fliplr(H_base));
H_fliplr(:, 2) = -H_fliplr(:, 2);
H = [H_base; H_fliplr];
inv_H = (H' * H) \ H';

%% - - - STBC - - - %%
stbc_temp = kron(reshape(info, 2, []), ones(1, 2));
stbc_temp(1 , 2 : 2 : end) = -conj(stbc_temp(1 , 2 : 2 : end));
stbc_temp(2 , 2 : 2 : end) = conj(stbc_temp(2 , 2 : 2 : end));

%% - - - ZF - - - %%
coded_info = stbc_temp;

%% - - - IQ Mod - - - %%
% CarrierSig_reshape = reshape(CarrierSig, sample_per_symbol, []);
% for counter = 1 : length(CarrierSig_reshape)
%    IQ_zf_mod_temp(1, :, counter) = coded_info(1, counter) * CarrierSig_reshape(:, counter);
%    IQ_zf_mod_temp(2, :, counter) = coded_info(2, counter) * CarrierSig_reshape(:, counter);
% end
% 
% IQ_zf_mod(1, :) = real(reshape(IQ_zf_mod_temp(1, :, :), 1, []));
% IQ_zf_mod(2, :) = real(reshape(IQ_zf_mod_temp(2, :, :), 1, []));
% 
% %% - - - IQ Demod - - - %%
% rx_data.one = reshape(IQ_zf_mod(1, :), 100, []);
% rx_data.two = reshape(IQ_zf_mod(2, :), 100, []);
% rx_zf_info.one = demodIQ(deCos, deSin, rx_data.one, sample_per_symbol);
% rx_zf_info.two = demodIQ(deCos, deSin, rx_data.two, sample_per_symbol);

%% - - - QAM Demod - - - %%
% rx_zf_qam_demod(1, :) = qamdemod(rx_zf_info.one, QAM_order);
% rx_zf_qam_demod(2, :) = qamdemod(rx_zf_info.two, QAM_order);

%% - - - Rx - - - %%
for counter_i = 1 : length(stbc_temp)
    rx_stbc(:, counter_i) = H_base * stbc_temp(:, counter_i);
end

%% - - -  STBC Estimate - - - %%
% for counter_i = 1 : length(stbc_temp)
%     demod_stbc(:, counter_i) = inv_H * rx_stbc(:, counter_i);
% end
rx_dt = rx_stbc(:, 1 : 2 : end);
rx_dt = [rx_dt; conj(rx_stbc(:, 2 : 2 : end))];

stbc_est = inv_H * rx_dt;
stbc_est(2, :) = conj(stbc_est(2, :));
data_de = zeros(1, 2 * length(stbc_est));
data_de(1: 2: end) = stbc_est(1, :);
data_de(2: 2: end) = stbc_est(2, :);

figure(1)
plot(real(rx_stbc), imag(rx_stbc), 'co');
hold on;
plot(real(data_de), imag(data_de), 'r*');
hold on;
plot(real(stbc_temp),imag(stbc_temp),'+')