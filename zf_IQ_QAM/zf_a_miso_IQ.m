clear all;
close all;
clc
% 2发1收的MISO，发射采用 Alamouti STBC ，zf预编码后调制输出
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
randi_bit = de2bi(randi_dec, k);
info = qammod(randi_dec, QAM_order);

%% - - - Channel- - - %%


%% - - - STBC - - - %%
stbc_temp = kron(reshape(info, 2, []), ones(1, 2));
stbc_temp(1 , 2 : 2 : end) = -conj(stbc_temp(1 , 2 : 2 : end));
stbc_temp(2 , 2 : 2 : end) = conj(stbc_temp(2 , 2 : 2 : end));

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

%% - - - IQ Demod - - - %%
Rx_data = reshape(IQ_zf_mod(1, :), 100, []);
rx_zf_info = demodIQ(deCos, deSin, Rx_data, sample_per_symbol);

%% - - - QAM Demod - - - %%
rx_zf_qam_demod = qamdemod(rx_zf_info, QAM_order);
rx_zf_bit = de2bi(rx_zf_qam_demod, k);

[numErrors_zf, ber_zf] = biterr(rx_zf_bit(:), randi_bit(:));
