clear all;
close all;
clc;

% 信道估计
%尝试使用LFM作为估计信号，在接收端做自相关获取信道
%的信道响应，根据冲击响应的值得到估计信道HEst
plot_flag = 0;
jay = sqrt(-1);
f1 = 10e3;
f2 = 20e3;
t = 0.02048;
fs = 50e3;
fraction =  64;
frequencyXaxis = (0 : pi/fraction : pi);
frequencyXaxis = frequencyXaxis(1 : 64)';
eta = 1e-3;    % Learning rate for LMS
H_base = [0.8, 0.5 + 0.1j, 0.3+ 0.2j, 0.6 + 0.4j];

pathDelays = [0 200 800 1200 2300 3700]*1e-9;    % sec
avgPathGains = [0 -0.9 -4.9 -8 -7.8 -23.9];      % dB
fD = 50;

%% - - - Gen LFM - - - %%
LFM_S = genLFM(f1, f2, t, fs);
% LFM_S = hilbert(LFM_S);

%% - - - Multipath - - - %%
% timeDelay(1, :) = floor(imag(H_base) * t * fs);   % zeropadding for delay
% timeDelay(2, :) = max(timeDelay(1, :)) - timeDelay(1, :); % zeropadding for align
% channelGain = real(H_base);
% rxLFMData = mulPathData(LFM_S, timeDelay, channelGain);

rayChan = comm.RayleighChannel(...
    'SampleRate',fs, ...
    'PathDelays',pathDelays, ...
    'AveragePathGains',avgPathGains, ...
    'NormalizePathGains',true, ...
    'RandomStream','mt19937ar with seed', ...
    'MaximumDopplerShift',fD, ...
    'Seed',22, ...
    'PathGainsOutputPort',true);
[rxLFMData, pathGains] = step(rayChan, LFM_S');
rxLFMData = rxLFMData';
%% - - - Do Something Special - - - %%

for counter_i = 1 : length(rxLFMData) / fraction
    LFM_SFraction(counter_i, :) = LFM_S((counter_i - 1) * fraction + 1 : counter_i * fraction);
    rxLFMDataFraction(counter_i, :) = rxLFMData((counter_i - 1) * fraction + 1 : counter_i * fraction);
end

fftLFM_S = fft(LFM_SFraction, [], 2);
fftrxLFM_S = fft(rxLFMDataFraction, [], 2);
%% LMS parameter
W = randn(1, fraction);                % Initial weights of LMS

    for n=1 : length(rxLFMData) / fraction

        y=W.*fftLFM_S(n,:);                       % Output of channel estimator
        e = fftrxLFM_S(n,:) - y;                   % Instantaneous error of LMS
        W = W + eta * e .* conj(fftLFM_S(n,:));   % Weight update rule of LMS
        J(n) = e * e';                          % Instantaneuous squared error
        
    end
    
    %% Results
% Impulse response
impulseResponse=invfreqz(conj(W'),frequencyXaxis,'complex',5,0);     % Estimated channel impulse response

Hest = LMSChannelEstimator(H_base);
%%
% 
% 
% % Wigner-Ville Distribution
% [ttfr,tt,tf] = wv(LFM_S);
% 
% [tfr,t,f] = wv(rxLFMData);
% 
% % Hough Transform
% [tht, trho, ttheta] = hough(ttfr, tf, tt);
% [ht, rho, theta] = hough(tfr, f, t);
% 
% % WHT Xcorr
%  result = xcorr2(tht, ht);

%% - - - Plot - - - %%
if plot_flag
    figure
    % f = f * (N1-1)/2/N1;
    t = t * 1/1000;
    [F, T] = meshgrid(f, t);
    mesh(F, T, abs(tfr));
    
    figure;
    mesh(theta*180/pi, rho, abs(ht));
    xlabel('theta'); ylabel('rho');
    
    figure;
    image(theta*(180/pi), rho, abs(ht));
    xlabel('theta'); ylabel('rho');
end
