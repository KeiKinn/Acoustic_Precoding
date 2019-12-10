clear all
close all
clc;
%% Instruction
% LS在时域估计，LMS在频域估计
% 本程序主要面向论文中第三章信道估计小节中对LS算法和LMS算法进行性能评估与对比

%% Basic Para
lineWidth = 2.2;
isPlot = 0;
MC = 10;
marker_color = [...
    0.0000    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840;...
    0.7500    0.7500    0.0000;...
    0.7500    0.0000    0.7500;...
    0.0000    0.5000    0.0000;...
    0.0000    0.0000    1.0000;...
    1.0000    0.0000    0.0000];

SNR = 6; %Signal to noise ratio
Rs = 63; %Let's define symbol rate for the plotting purposes
errorResult = [];
mc_errortemp = zeros(2, 31);
%% Channel creation and channel modelling
p = [0.19+.56j .45-1.28j -.14-.53j -.19+.23j .33+.51j]; % Example complex channel model
[H,f]=freqz(p,1,-Rs/2:1:Rs/2,Rs);
L1=1; % Channel maximum tap is the second one

%% Creation of the data
N = (2^10);    % Number of samples
Bits = 2;      % For modulation
for mc = 1 : MC
for SNR = -20:10
    index = SNR + 21;
    data = randi([0 ((2^Bits)-1)], 1, N);  % Random signal
    Ak = pskmod(data, 2^Bits);          % Phase shit keying (PSK) modulation
    % Ak =  lteZadoffChuSeq(1, N - 1)'; % ZC sequnence
    
    %% Signal Processing
    Rk = filter(p,1,Ak); % Received signal
    Rk = Rk(L1+1:end); % Discard the delay, if channel has pre-cursor taps
    noise = (1/sqrt(2))*(randn(size(Rk)) + 1j*randn(size(Rk))); %Initial noise vector
    P_s = var(Rk); % Signal power
    P_n = var(noise); % Noise power
    % Defining noise scaling factor based on the desired SNR:
    noise_scaling_factor = sqrt(P_s/P_n./10.^(SNR./10));
    Rk_noisy=Rk+noise*noise_scaling_factor; % Received signal
    Nsig = mdltest(Rk_noisy')
    %% LS
    M = 30; %number of reference symbols used for channel estimation; this is pretty small amount
    estimate_length = 7; %this defines how long is the channel estimate's impulse response
    A_conv=convmtx(Ak(2:M+1).',estimate_length); % Convolution matrix
    p_LS=((A_conv'*A_conv)\(A_conv'))*Rk_noisy(1:size(A_conv,1)).'; % LS solution
    %%% square error
    %     LSk = filter(p_LS, 1, Ak);
    %     e_LS = LSk(L1+1 : end) - Rk;
    e_LS = HE - H;
    errorResult(1, index) = std(abs(e_LS)) ;
    %% LMS
    eta =2.9e-4; % step size
    fft_Ak = fft(Ak, N);
    fft_Rk_noisy = fft(Rk, N);
    fft_Ak_frac = reshape(fft_Ak, 16, []);
    sys_out = fft_Ak_frac .* repmat(H, 16, 1);
    fft_Rk_noisy_frac = reshape(fft_Rk_noisy, 16, []);
    fft_Rk_noisy_frac = awgn(sys_out, SNR);
    
    W = ones(size(H));
    for n = 1 : 16
        y = W .* fft_Ak_frac(n, :);
        e = fft_Rk_noisy_frac(n, :) - y;
        W = W + eta * conj(fft_Ak_frac(n, :)) .*e;
        J(n) = e * e';
    end
    e_LMS = W - H;
    errorResult(2, index) =std(abs(e_LMS)) ;
    %% Plotting amplitude response of the channel:
    if isPlot
        figure(2)
        plot(20*log10(abs(H)), 'color', marker_color(1, :), 'LineWidth', lineWidth);
        hold on
        plot(20*log10(abs(HE)), 'color',marker_color(2, :), 'LineWidth', lineWidth);
        hold on
        plot(20*log10(abs(W)), 'color',marker_color(5, :), 'LineWidth', lineWidth);
        legend('Channel','LS Channel Estimate', 'LMS Channel Estimation');
        xlabel('Frequency[kHz] ');ylabel('Amplitude response [dB]');
        title('Campare Between LMS and LS on Channel Estimation in Channel Frequency Response');
        grid on
        figloc
    end
    
end
mc_errortemp = mc_errortemp + errorResult;
end
%% Error
figure('Name', 'error')
semilogy([-20:10], mc_errortemp(1,:)/10, 'LineWidth', lineWidth);
hold on
semilogy([-20:10], mc_errortemp(2,:)/10, 'LineWidth', lineWidth);
xlabel('SNR[dB]');ylabel('MSE');
legend('LS','LMS');
title('Compare MSE between LS and LSM')
grid on
figloc

%% end
% figure(3)
% % subplot(211)
% stem(-L1:length(p)-L1-1,abs(p),'b', 'LineWidth', lineWidth);
% grid on;
% hold on;
% % LS result
% stem(-L1:length(p_LS)-L1-1,abs(p_LS),'r', 'LineWidth', lineWidth);
% legend('Channel','LS channel estimate');
% title('Absolute values of the impulse responses')
% hold on;
% % LMS result
% stem(-L1:length(P_LMS)-L1-1,abs(P_LMS),'k', 'LineWidth', lineWidth);
% legend('Channel','LS channel estimate', 'LMS channel estimate');
% title('Absolute values of the impulse responses')
% % subplot(212)
% % plot(20*log10(abs(H)), 'color', marker_color(1, :), 'LineWidth', lineWidth);
% % xlabel('Frequency[kHz] ');ylabel('Amplitude response [dB]');
% % title('Values of Frequency Response')
% % grid on
% % hold off;
% figloc

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

