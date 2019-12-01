% 本程序主要面向论文中第三章预编码方案中的mmse和zf预编码方案的性能对比
clear all
close all
clc;

format long; 
Nt=2;
Nr=2;
SNR=[-20:2:16];
channel_n=100*ones(1,length(SNR));
error_mmselinp=zeros(1,length(SNR));
error_zflinp=zeros(1,length(SNR));
for loop_ebno=1:length(SNR)
    snr=10.^(SNR(loop_ebno)/10);
    ea=1;
    es=ea*Nt;
    sigma_n2=es/snr;
    num=200;
    for loop_channel=1:channel_n(loop_ebno)
        H=sqrt(1/2)*(randn(Nr,Nt)+j*randn(Nr,Nt));
        mmse_F=H'/(H*H'+Nt/snr*eye(Nt));
        zf_F = H'/(H*H');
        beta_mmse=sqrt(es/norm(mmse_F,'fro').^2);
        beta_zf=sqrt(es/norm(zf_F,'fro').^2);
        F_mmse=beta_mmse*mmse_F; 
        F_zf=beta_zf*zf_F;
        for loop_num=1:num
            gen_u=(sign(randn(Nt,1))+j*sign(randn(Nt,1)));
            u=sqrt(1/2)*gen_u;
            x_mmse=F_mmse*u;
            x_zf=F_zf*u;
            noise=sqrt(sigma_n2/2)*(randn(Nr,1)+j*randn(Nr,1));
            noise1=sqrt(sigma_n2/2)*(randn(Nr,1)+j*randn(Nr,1));
            y_mmse=H*x_mmse+noise;
            y_zf=H*x_zf+noise1;
            r_mmse=1/beta_mmse*y_mmse;
            r_zf=1/beta_zf*y_zf;
            rev_data_mmse=sign(real(r_mmse))+j*sign(imag(r_mmse));
            rev_data_zf=sign(real(r_zf))+j*sign(imag(r_zf));
            error_mmselinp(1,loop_ebno)=error_mmselinp(1,loop_ebno)+sum(((abs(rev_data_mmse-gen_u)).^2)/4);
            error_zflinp(1,loop_ebno)=error_zflinp(1,loop_ebno)+sum(((abs(rev_data_zf-gen_u)).^2)/4);
        end 
    end 
    ber_mmselinp(1,loop_ebno)=error_mmselinp(1,loop_ebno)/(num*Nt*2*channel_n(loop_ebno));
                                                                                          
    ber_zflinp(1,loop_ebno)=error_zflinp(1,loop_ebno)/(num*Nt*2*channel_n(loop_ebno));
end

P1=semilogy(SNR,ber_mmselinp);
hold on
P2=semilogy(SNR,ber_zflinp);
set(P1,'Linewidth',[2]);
set(P2,'Linewidth',[2]);
grid on;
xlabel('Symbol SNR(dB)');ylabel('BER');
title('Compare between MMSE and ZF precoder')
leg1='MMSE Precoder';
leg2='ZF Precoder';
legend(leg1,leg2);
figloc;


