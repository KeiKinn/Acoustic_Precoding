clear all;
close all;
clc;

% 信道估计
%尝试使用LFM作为估计信号，在接收端做自相关获取信道
%的信道响应，根据冲击响应的值得到估计信道HEst
plot_flag = 1;
jay = sqrt(-1);
f1 = 10e3;
f2 = 20e3;
t = 0.02;
fs = 50e3;
H_base = [0.8, 0.5 + 0.1j, 0.3+ 0.2j, 0.6 + 0.4j];

pathDelays = [0 200 800 1200 2300 3700]*1e-9;    % sec
avgPathGains = [0 -0.9 -4.9 -8 -7.8 -23.9];      % dB
fD = 50;

%% - - - Gen LFM - - - %%
LFM_S = genLFM(f1, f2, t, fs);
LFM_S = hilbert(LFM_S);

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

rxLFMData = rayChan(LFM_S');
rxLFMData = rxLFMData';
%% - - - Do Something Special - - - %%
% Wigner-Ville Distribution
[ttfr,tt,tf] = wv(LFM_S);

[tfr,t,f] = wv(rxLFMData);

% Hough Transform
[tht, trho, ttheta] = hough(ttfr, tf, tt);
[ht, rho, theta] = hough(tfr, f, t);

% WHT Xcorr
 result = xcorr2(tht, ht);

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
