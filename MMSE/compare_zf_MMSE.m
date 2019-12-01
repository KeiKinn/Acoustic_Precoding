% 本程序主要面向论文中第三章预编码方案中的mmse和zf预编码方案的性能对比
% 2019.11.30 尝试重构程序，主要关注代码格式，信号形式，信道，性能评估以及性能对比
clear all
close all
clc;

%% Basic Para
isPlot = 0;
MCNum = 1;
Nt=2; % transmitter number
Nr=2; % receiver number
ea=1;
es=ea*Nt;
SNR=[-20:0.1:0];
% SNR = -4;
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
    for n = 1 : length(SNR)
        %         %% Channel
        %         H=sqrt(1/2)*(randn(Nr,Nt)+1i*randn(Nr,Nt));
        
        %% MMSE Precoding
        mmse_F=H'/(H*H'+Nt/snr(n)*eye(Nt));
        beta_mmse=sqrt(es/norm(mmse_F,'fro').^2);
        F_mmse=beta_mmse*mmse_F;
        %% ZF Precoding
        zf_F = H'/(H*H');
        beta_zf=sqrt(es/norm(zf_F,'fro').^2);
        F_zf=beta_zf*zf_F;
        
        %% Signal in Transmitter
        [info, randi_bit] = genQAM(qamOrder, symbolNum);
        STBCData = genSTBC(info);
        signalMMSE = F_mmse * STBCData;
        signalZF = F_zf * STBCData;
        
        %% Through Channel
        dataSTBCTemp = H * STBCData;    % without precoding
        dataSTBC = awgn(dataSTBCTemp, SNR(n));
        
        dataMMSETemp = H * signalMMSE;  % Channel State
        dataMMSE = awgn(dataMMSETemp, SNR(n)); % noise
        
        dataZFTemp = H * signalZF;
        dataZF = awgn(dataZFTemp, SNR(n));
        
        %% Receiver
        H_pseSTBC = H;
        invPseSTBC = invPseH(H_pseSTBC);    % Psedo Inverse for STBC   
        rxSTBC = shapeMat(dataSTBC);
        deSTBC = invPseSTBC * rxSTBC;
        
        H_pseMMSE = H * F_mmse;
        invPseMMSE = invPseH(H_pseMMSE); % Psedo Inverse for MMSE-STBC   
        rxMMSE = shapeMat(dataMMSE);
        deMMSE = invPseMMSE * rxMMSE; % Decode STBC
        
        H_pseZF = H * F_zf;
        invPseZF = invPseH(H_pseZF);      
        rxZF = shapeMat(dataZF);
        deZF = invPseZF * rxZF;
        %% BER
        de_bit_STBC = qam2bit(deSTBC, qamOrder);
        de_bit_MMSE= qam2bit(deMMSE, qamOrder);
        de_bit_ZF = qam2bit(deZF, qamOrder);
        
        [numErrorsSTBC, berSTBC(n)] = biterr(de_bit_STBC(:), randi_bit(:));
        [numErrorsMMSE, berMMSE(n)] = biterr(de_bit_MMSE(:), randi_bit(:));
        [numErrorsZF, berZF(n)] = biterr(de_bit_ZF(:), randi_bit(:));
    end
    
    tempSTBC = tempSTBC + berSTBC;
    tempMMSE = tempMMSE + berMMSE;
    tempZF = tempZF + berZF;
end
%% end

if isPlot
    figure('Name', 'Constellation: SNR = -5dB')
    plot(real(deMMSE(:)), imag(deMMSE(:)), 'co');
    hold on;
    plot(real(deZF(:)), imag(deZF(:)), 'r*');
    hold on;
    plot(real(info), imag(info),'ko', 'MarkerSize',6 , 'MarkerFaceColor','k');
    hold off;
    legend("MMSE Precoder", "ZF Precoder", "16 QAM", 'FontSize', 16);
    xlabel('RE'); ylabel('IM');
    title('Constellation Compare bwtween MMSE Precoder and ZF Precoder with SNR = -5dB');
    figloc
    
    figure('Name', 'MMSE Constellation: SNR = -5dB')
    plot(real(deMMSE(:)), imag(deMMSE(:)), 'co');
    hold on;
    plot(real(info), imag(info),'ko', 'MarkerSize',6 , 'MarkerFaceColor','k');
    hold off;
    legend("MMSE Precoder", "16 QAM", 'FontSize', 16);
    xlabel('RE'); ylabel('IM');
    title('Constellation of MMSE Precoder with SNR = -5dB');
    figloc
    
    figure('Name', 'ZF Constellation: SNR = -5dB')
    plot(real(deZF(:)), imag(deZF(:)), 'r*');
    hold on;
    plot(real(info), imag(info),'ko', 'MarkerSize',6 , 'MarkerFaceColor','k');
    hold off;
    legend("ZF Precoder", "16 QAM", 'FontSize', 16);
    xlabel('RE'); ylabel('IM');
    title('Constellation of ZF Precoder with SNR = -5dB');
    figloc
    
    figure('Name', 'MC ber')
    semilogy(SNR, berSTBC / MCNum);
    hold on;
    semilogy(SNR, berMMSE / MCNum);
    hold on;
    semilogy(SNR, berZF / MCNum);
    hold off;
    grid on;
    xlabel('Symbol SNR(dB)');ylabel('BER');
    legend('Only STBC', 'MMSE Precoder','ZF Precoder')
    title('Compare among only STBC, MMSE and ZF precoder')
    figloc;
end