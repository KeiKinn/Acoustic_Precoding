clear all;
close all;
clc
% 将数字信号通过IQ调制的方式加载到载波上，并通过多径信道，未叠加噪声
fc = 10e3;
fs = 50e3;
QAM_protype = repmat(-3:2:3,1,4)-[repmat(-3j,1,4) repmat(-1j,1,4) repmat(1j,1,4) repmat(3j,1,4)];
% info = QAM_protype(randi(16,1000,1)); % 16-QAM sequence of 1000 samples
% 
sample_per_time = 1000;
jay = sqrt(-1);
deltaT = 1 / fs;
sample_per_symbol = 100;
info = QAM_protype(randi(16,sample_per_time,1)); % 16-QAM sequence of ,sample_per_time samples
N1 = 0; N2=8;   % equalizer sum_(-N1)^(N2)
L1 = 0; L2 = 4;  % isi sum_(-L1)^(N2)
info = [info zeros(1, L1 + N1)];

channel_gain = [0.69, 0.18, 0.05, 0.06, 0.02];
channel_delay = [0, 16, 22, 28, 41] * 1e-4;

symbol_rate = fs / sample_per_symbol;
SigTim = sample_per_time / symbol_rate;
DataTim = 0 : deltaT : SigTim - deltaT;
CarrierSig = exp(jay * 2 * pi * fc * DataTim);
deCos = cos(2 * pi * fc * (0 : sample_per_symbol - 1) / fs);
deSin =  -sin(2 * pi * fc * (0 : sample_per_symbol - 1) / fs);

actual_isi = channel_gain + jay .* channel_delay;
length_equalizer = N1 + N2 + 1; 
length_actual_isi = length(actual_isi);
time_delay(1, :) = floor(imag(actual_isi) * fs);   % zeropadding for delay
time_delay(2, :) = max(time_delay(1, :)) - time_delay(1, :); % zeropadding for align

%%- - - ZF Precoding - - -%%
destina_matrix = zeros(L1 + L2 + N1 + N2 + 1, 1);
destina_matrix(L1 + N1 + 1) = 1;
matrix_isi = convmtx(actual_isi.', length_equalizer);
matrix_zf = ((matrix_isi' * matrix_isi) \ (matrix_isi')) * destina_matrix; % (P^HP)^(-1)P^H

coded_info = filter(matrix_zf, 1, info);
coded_info = coded_info(L1 + N1 + 1 : end);

%%- - - IQ Modulation - - -%%
CarrierSig_reshape = reshape(CarrierSig, sample_per_symbol, []);
for counter = 1 : length(CarrierSig_reshape)
   IQ_info_modu(:, counter) = info(counter) * CarrierSig_reshape(:, counter);
   IQ_zf_modu(:, counter) = coded_info(counter) * CarrierSig_reshape(:, counter);
end
IQ_info_modu = real(reshape(IQ_info_modu, 1, []));
IQ_zf_modu = real(reshape(IQ_zf_modu, 1, []));

%%- - - Muliti Path - - -%%
%%%- - - info - - -%%%
rx_info_data = mulPathData(IQ_info_modu, time_delay, channel_gain);
rx_info_data = rx_info_data(1 : length(DataTim));
Rx_info_data = reshape(rx_info_data, 100, []);
%%%- - - zf_info - - -%%%
rx_data = mulPathData(IQ_zf_modu, time_delay, channel_gain);
rx_data = rx_data(1 : length(DataTim));
Rx_data = reshape(rx_data, 100, []);

%%- - - demod - - -%%
rx_info = demodIQ(deCos, deSin, Rx_info_data, sample_per_symbol);
rx_zf_info = demodIQ(deCos, deSin, Rx_data, sample_per_symbol);

%%- - - Plot - - -%%
figure(1)
plot(real(rx_info), imag(rx_info), 'cd')
hold on;
plot(real(coded_info), imag(coded_info), 'r*')
hold on;
plot(real(rx_zf_info), imag(rx_zf_info), 'o');
hold on;
plot(real(info), imag(info),'ko', 'MarkerSize',6 , 'MarkerFaceColor','k');
grid on;
legend("ISI Channel without ZF", "ZF Coded","ISI Channel with ZF", "Ideal constellation points");
xlabel('RE'); ylabel('IM');
figloc

% figure(2)
% subplot(2, 1, 1)
% stem(-L1 : L2, abs(actual_isi), 'b', 'LineWidth', 1.3);
% legend('ISI Channel');
% grid on;
% title('Absolute values of impulse responses'); % Absolute values of channel impulse response
% subplot(2,1,2)
% stem(-N1-L1 : N2+L2, abs(conv(actual_isi, matrix_zf)), 'LineWidth', 1.3); 
% legend('ISIS Channel + ZF Precoding');
% grid on;
% title('Absolute values of impulse responses'); % Absolute values of channel impulse response
