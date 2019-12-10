% 本程序主要面向论文中第三章预编码方案中的mmse和zf预编码方案的性能对比
% 2019.11.30 尝试重构程序，主要关注代码格式，信号形式，信道，性能评估以及性能对比
clear all
close all
clc;

%% Basic Para
isPlot = 0;
lineWidth = 2.2;
marker_style = {'o-','s--','v-.','+:','<-','>--','x-.','^:','*-','d--','h-.','p:'};
MCNum = 1;
Nt=1; % transmitter number
Nr=1; % receiver number
ea=1;
es=ea*Nt;
% SNR=[-20:2:0];
SNR = -6;
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
H = [0.16+0.34j, 0.65-1.48j, -0.14-.93j, .43+.23j, 0.453+.51j;
        0.19+0.56j, 0.45-1.28j, -0.14-.53j, -0.19+.23j, 0.33+.51;
        0.13-0.45j, 0.67+0.58j, 0.44+.43j, -0.14-.13j, 0.153+.91];
for MC = 1 : MCNum
    if ~mod(MC, 100)
        display(floor(MC/100));
    end
    for n = 1 : length(SNR)
        %% ZF Precoding
        zf_F = H'/(H*H');
        beta_zf=sqrt(es/norm(zf_F,'fro').^2);
        F_zf=beta_zf * zf_F;
        
        %% Signal in Transmitter
        [info, randi_bit] = genQAM(qamOrder, symbolNum);
        infoTemp = [info info info]';
        
        signalZF = F_zf * infoTemp;
        
        %% Through Channel
        
        
        dataZFTemp = H * signalZF;
        dataZF = awgn(dataZFTemp, SNR(n));
        
        %% Receiver
        
        
        H_pseZF = H * F_zf;
        invPseZF = invPseH(H_pseZF);
        rxZF = shapeMat(dataZF);
        deZF = invPseZF * rxZF;
        
        %% BER
        
        de_bit_ZF = qam2bit(deZF, qamOrder);
        
        
        [numErrorsZF, berZF(n)] = biterr(de_bit_ZF(:), randi_bit(:));
    end
    
    tempZF = tempZF + berZF;
end
%% end

if isPlot
    
end