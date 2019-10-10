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
H_base = [0.8, 0.5 + 0.1j, 0.3 + 0.2j, 0.6 + 0.4j];

%% - - - Gen LFM - - - %%
LFM_S = genLFM(f1, f2, t, fs);

%% - - - Multipath - - - %%
timeDelay(1, :) = floor(imag(H_base) * t * fs);   % zeropadding for delay
timeDelay(2, :) = max(timeDelay(1, :)) - timeDelay(1, :); % zeropadding for align
channelGain = real(H_base);
rxLFMData = mulPathData(LFM_S, timeDelay, channelGain);

%% - - - Do Something Special - - - %%
% Wigner-Ville Distribution
[ttfr,tt,tf] = wv(LFM_S);

[tfr,t,f] = wv(rxLFMData);

% Hough Transform
[tht, trho, ttheta] = hough(ttfr, tf, tt);
[ht, rho, theta] = hough(tfr, f, t);

% WHT Xcorr
% result = xcorr2(tht, ht);

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
    
end
