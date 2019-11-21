clear all;
close all;
clc
%产生N点 16QAM 信号，经过 ZF 后通过 ISI 信道，通过查看星座图估计算法效果

N = 1000;
QAM_protype = repmat(-3:2:3,1,4)-[repmat(-3j,1,4) repmat(-1j,1,4) repmat(1j,1,4) repmat(3j,1,4)];
info = QAM_protype(randi(16,N,1)); % 16-QAM sequence of ,sample_per_time samples
N1 = 3; N2=5;   % equalizer sum_(-N1)^(N2)
L1 = 1; L2 = 3;  % isi sum_(-L1)^(N2)
SNR = 5;
snr=10.^(SNR/10);
%%- - - Channel - - -%%
actual_isi = [0.19+0.56j, 0.45-1.28j, -0.14-.53j, -0.19+.23j, 0.33+.51j];
length_equalizer = N1 + N2 + 1; 
length_actual_isi = length(actual_isi);

%%- - - ZF Precoding - - -%%
matrix_isi = convmtx(actual_isi.', length_equalizer);
destina_matrix = zeros(L1 + L2 + N1 + N2 + 1, 1);
destina_matrix(L1 + N1 + 1) = 1;
matrix_zf = ((matrix_isi' * matrix_isi) \ (matrix_isi')) * destina_matrix; % (P^HP)^(-1)P^H pesuoinverse
[~, c] = size(matrix_isi);
matrix_zf = ((matrix_isi' * matrix_isi + eye(c)/snr) \ (matrix_isi')) * destina_matrix;
coded_info = filter(matrix_zf, 1, info);
coded_info = coded_info(L1 + N1 + 1 : end);

%%- - - Rx_Signal - - -%%
noise = sqrt(0.5/snr) * (randn(N, 1) + 1i * randn(N, 1));
rxPreSig = coded_info + noise(L1 + N1 + 1 : end)';
rxSig = info + noise';
Rx_info = filter(actual_isi, 1, rxSig);
Rx_zf_info = filter(actual_isi, 1, rxPreSig);

%%- - - Plot - - -%%
figure
plot(real(Rx_info), imag(Rx_info), 'o');
hold on;
plot(real(Rx_zf_info), imag(Rx_zf_info), 'o');
hold on;
plot(real(info), imag(info),'ko', 'MarkerSize',5 , 'MarkerFaceColor','k');
grid on;
legend("ISI Channel without MMSE", "ISI Channel with MMSE", "Ideal constellation points");
figloc;