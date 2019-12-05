% 本程序主要面向论文中第三章预编码方案中的mmse和zf预编码方案的性能对比
% 2019.11.30 尝试重构程序，主要关注代码格式，信号形式，信道，性能评估以及性能对比
clear all
close all
clc;

%% Basic Para
isPlot = 0;
qamSymbol = [...
    -3-3i,-3-1i,-3+3i,-3+1i, ...
    -1-3i,-1-1i,-1+3i,-1+1i, ...
    +3-3i,+3-1i,+3+3i,+3+1i, ...
    +1-3i,+1-1i,+1+3i,+1+1i ];
lineWidth = 2.2;
MCNum = 1;
Nt=2; % transmitter number
Nr=2; % receiver number
ea=1;
es=ea*Nt;
% SNR=[-20:2:0];
SNR = -4;
snr = 10.^(SNR/10); % sigma_s / sigma_n for MMSE
qamOrder = 16;
symbolNum = 1e3;
channel_n=100*ones(1,length(SNR));

%% Initial Martix
berSTBC = zeros(1, length(SNR));
berMMSE=zeros(1,length(SNR));
berZF=zeros(1,length(SNR));
tempSTBC = zeros(1, length(SNR));
tempMMSE = zeros(1, length(SNR));
tempZF =zeros(1, length(SNR));
%% Channel
H=sqrt(1/2)*(randn(Nr,Nt)+1i*randn(Nr,Nt));

for MC = 1 : MCNum
    %% Counter
    if ~mod(MC, 100)
        clc;
        display(floor(MC/100));
    end
    
    for n = 1 : length(SNR)
        %% MMSE Precoding
        mmse_F = H'/(H*H'+Nt/snr(n)*eye(Nt));
        beta_mmse = sqrt(es/norm(mmse_F,'fro').^2);
        F_mmse = beta_mmse*mmse_F;
        
        %% ZF Precoding
        zf_F = H'/(H*H');
        beta_zf=sqrt(es/norm(zf_F,'fro').^2);
        F_zf=beta_zf * zf_F;
        
        %% Signal in Transmitter
        [info, randi_bit] = genQAM(qamOrder, symbolNum);
        STBCData = genSTBC(info);
        STBCData = reshape(info, 2, []);
        
        for counter = 1 : symbolNum/2
            [signalSQUID(:, counter), beta_squid(counter)] = squid(STBCData(:, counter), H, snr);
        end
        beta_squid = [beta_squid; beta_squid];
        %% Through Channel

        dataSQUIDTemp = H * signalSQUID;
        dataSQUID = beta_squid .* awgn(dataSQUIDTemp, SNR);
        
        %% Receiver
% %         invPseSQUID = eye(2, 4);
% %         rxSQUID = shapeMat(dataSQUID);
% %         deSQUID = invPseSQUID * rxSQUID;
deSQUID = dataSQUID;
        %% BER

        de_bit_SQUID = qam2bit(reshape(deSQUID, 2, []), qamOrder);
        
        [numErrorsSQUID, berSQUID(n)] = biterr(de_bit_SQUID(:), randi_bit(:));
        
    end
end
%% end

if isPlot
    
end