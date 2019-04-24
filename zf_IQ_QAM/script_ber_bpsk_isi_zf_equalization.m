% Script for computing the BER for BPSK modulation in 3 tap ISI 
% channel. Zero Forcing equalization with 3/5/7/9 tap is performed
% and the BER computed

clear
N  = 10^6; % number of bits or symbols

Eb_N0_dB = [0:10]; % multiple Eb/N0 values
K = 4;


for ii = 1:length(Eb_N0_dB)

   % Transmitter
   ip = rand(1,N)>0.5; % generating 0,1 with equal probability
   s = 2*ip-1; % BPSK modulation 0 -> -1; 1 -> 0 

   % Channel model, multipath channel
   nTap = 3;
   ht = [0.2 0.9 0.3]; 

   chanOut = conv(s,ht);  
   n = 1/sqrt(2)*[randn(1,N+length(ht)-1) + j*randn(1,N+length(ht)-1)]; % white gaussian noise, 0dB variance 
   
   % Noise addition
   y = chanOut + 10^(-Eb_N0_dB(ii)/20)*n; % additive white gaussian noise

   for kk = 1:K

     L  = length(ht);
     hM = toeplitz([ht([2:end]) zeros(1,2*kk+1-L+1)], [ ht([2:-1:1]) zeros(1,2*kk+1-L+1) ]);
     d  = zeros(1,2*kk+1);
     d(kk+1) = 1;
     c  = [inv(hM)*d.'].';

     % mathched filter
     yFilt = conv(y,c);
     yFilt = yFilt(kk+2:end); 
     yFilt = conv(yFilt,ones(1,1)); % convolution
     ySamp = yFilt(1:1:N);  % sampling at time T
   

     % receiver - hard decision decoding
     ipHat = real(ySamp)>0;

     % counting the errors
     nErr(kk,ii) = size(find([ip- ipHat]),2);

   end
   

end

simBer = nErr/N; % simulated ber
theoryBer = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10))); % theoretical ber

% plot
close all
figure
semilogy(Eb_N0_dB,simBer(1,:),'bs-'),'Linewidth',2;
hold on
semilogy(Eb_N0_dB,simBer(2,:),'gd-'),'Linewidth',2;
semilogy(Eb_N0_dB,simBer(3,:),'ks-'),'Linewidth',2;
semilogy(Eb_N0_dB,simBer(4,:),'mx-'),'Linewidth',2;
axis([0 10 10^-3 0.5])
grid on
legend('sim-3tap', 'sim-5tap','sim-7tap','sim-9tap');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for BPSK in ISI with ZF equalizer');
figloc
