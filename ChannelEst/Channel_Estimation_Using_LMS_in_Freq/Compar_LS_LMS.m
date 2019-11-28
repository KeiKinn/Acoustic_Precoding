clear all
close all
clc;

SNR = 15; %Signal to noise ratio
Rs = 63; %Let's define symbol rate for the plotting purposes

%% Creation of the data

N = (2^10);    % Number of samples
Bits = 2;      % For modulation
data = randi([0 ((2^Bits)-1)], 1, N);  % Random signal
Ak = pskmod(data, 2^Bits);          % Phase shit keying (PSK) modulation
% Ak =  lteZadoffChuSeq(1, N - 1)'; % ZC sequnence

%% Channel creation and channel modelling
p = [0.19+.56j .45-1.28j -.14-.53j -.19+.23j .33+.51j]; % Example complex channel model
[H,f]=freqz(p,1,-Rs/2:1:Rs/2,Rs);
L1=1; % Channel maximum tap is the second one

Rk = filter(p,1,Ak); % Received signal
Rk = Rk(L1+1:end); % Discard the delay, if channel has pre-cursor taps

noise = (1/sqrt(2))*(randn(size(Rk)) + 1j*randn(size(Rk))); %Initial noise vector
P_s = var(Rk); % Signal power
P_n = var(noise); % Noise power
% Defining noise scaling factor based on the desired SNR:
noise_scaling_factor = sqrt(P_s/P_n./10.^(SNR./10));
Rk_noisy=Rk+noise*noise_scaling_factor; % Received signal

M = 30; %number of reference symbols used for channel estimation; this is pretty small amount
estimate_length = 7; %this defines how long is the channel estimate's impulse response
A_conv=convmtx(Ak(2:M+1).',estimate_length); % Convolution matrix
p_LS=((A_conv'*A_conv)\(A_conv'))*Rk_noisy(1:size(A_conv,1)).'; % LS solution

% Plotting amplitude response of the channel:
figure(2)
plot(20*log10(abs(H)),'b');
hold on
figure(2)
[HE,fE]=freqz(p_LS,1,-Rs/2:1:Rs/2,Rs);
plot(20*log10(abs(HE)),'r'); legend('Channel','LS Channel Estimate');
xlabel('Frequency ');ylabel('Amplitude response [dB]'); legend('Channel');
figloc

% figure(3)
% stem(-L1:length(p)-L1-1,abs(p),'k');
% hold on;
% stem(-L1:length(p_LS)-L1-1,abs(p_LS),'r');
% legend('Channel','LS channel estimate');
% title('Absolute values of the impulse responses')
% hold off;
% figloc

%% LMS
eta =2.7e-4;

fft_Ak = fft(Ak, N);
fft_Rk_noisy = fft(Rk, N);
fft_Ak_frac = reshape(fft_Ak, 16, []);
fft_Rk_noisy_frac = reshape(fft_Rk_noisy, 16, []);

W = randn(size(H));
for n = 1 : 16
    y = W .* fft_Ak_frac(n, :);
    e = fft_Rk_noisy_frac(n, :) - y;
    W = W + eta * e .* conj(fft_Ak_frac(n, :));
    J(n) = e * e';
end
figure(2)
% plot(J);
% H_Est=invfreqz(conj(W'), f, 'complex', 4, 0);
plot(20*log10(abs(flip(W))), 'g');

% %% MSE Equalizer (example with 10000 training symbols)
% % Initialization
% FII = zeros(31,31); % Autocorrelation Matrix initialization
% alfa = zeros(31,1); % Cross-correlation vector initialization
% % Estimating FII and alfa using sample estimates based on training symbols
% % Notice that here we use all the generated data as training symbols
% for i = 16:length(Rk_noisy)-15,
%     rk = flipud(Rk_noisy(i-15:i+15).'); % Received signal vector
%     FII = FII + rk*conj(rk).'; % Autocorrelation matrix
%     alfa = alfa + Ak(i)*conj(rk);
% end
% FII = FII/(length(Rk_noisy)-30); % Final sample estimate of the autocorrelation matrix
% alfa = alfa/(length(Rk_noisy)-30); % Final sample estimate of the cross-correlation vector
% c_MSE = inv(conj(FII))*alfa; % Equalizer coefficients
% 
% figure(4); stem(abs(conv(p,c_MSE)));
% title('Effective impulse reponse (abs) of the MSE equalized system')
% figloc
% 
% figure(5); [H,f]=freqz(p,1,-Rs/2:Rs/200:Rs/2,Rs); plot(f/1e6,20*log10(abs(H)),'b', 'LineWidth', lineWidth);
% xlabel('Frequency');ylabel('Amplitude response');
% figure(5); hold on; [Hf,f]=freqz(c_MSE,1,-Rs/2:Rs/200:Rs/2,Rs);
% plot(f/1e6,20*log10(abs(Hf)),'k','LineWidth', lineWidth);
% hold on;
% H_t = H.*Hf;
% plot(f/1e6,20*log10(abs(H_t)),'r', 'LineWidth', lineWidth); hold off;
% legend('ISI Channel','MMSE Precoder', 'Channer after MMSE');
% grid on;
% figloc

