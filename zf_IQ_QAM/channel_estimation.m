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

%% - - - Plot - - - %%
if plot_flag
    fftp(rxLFMData, fs);
end
