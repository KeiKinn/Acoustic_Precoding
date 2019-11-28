%% Channel Estimation Using Least Mean Square (LMS) algorithm in Frequency domain
% Author: SHUJAAT KHAN

% Channel estimation or system identification is a technique in which we model 
% the relationship between given inputs and outputs of the unknown system.

%% Start
clc;
clear all;
close all;

%% Initializing simulation parameters
runs = 1;     % Number of Monte Carlo simulations (runs)
temp1 = 0;     % temporary estimated coefficients of channel (for each run)
temp2 = 0;     % tempporary mean squared error (for each run) 
N = (2^10);    % Number of samples
Bits = 2;      % For modulation    
fft_len=64;    % Fast Fourier Transform (Frame size) or channel length
SNR = 10;      % Signal to noise ratio or noise level
eta = 1e-2;    % Learning rate for LMS

%% Defining Unknown channel
% Channel impulse response
channel_impulse_response = [0.1+0.2j 0.2 0.5 -0.7];    
% Frequency response of channel
[com_freq,Freq_vector]=freqz(channel_impulse_response,1,fft_len);

%% Monte Carlos Simulation
for k = 1 : runs
    % Generate Random signal for each independent run
    data = randi([0 ((2^Bits)-1)],1,N);  % Random signal
    x = pskmod(data,2^Bits);          % Phase shit keying (PSK) modulation
    % x = qammod(data,2^Bits);        % Quadrature amplitude modulation (QAM)

    % zero padding to ensure each frame is of size equals to the channel/frame
    padding=rem(length(x),fft_len);
    if padding~=0
    padding=fft_len-padding;
    end
    x=[x zeros(1,padding)];
    signal_length=length(x); % final length after padding

    % Dividing data stream into frames of length equals to fft_len
    x=reshape(x,[signal_length/fft_len fft_len]);   
    
    % Convert time domain signal to frequency domain signal using FFT
    signal=fft(x,[],1);
    
    % Generate desired output signal by (Convolution in time domain or
    % multiplication in frequency domain)
    sys_out=signal.*repmat(conj(com_freq'),signal_length/fft_len,1);    % signal after passing through channel

desired = sys_out;  % Addition of white gaussian noise (desired signal after disturbance)

%% LMS parameter
W = randn(size(com_freq'));                % Initial weights of LMS

    for n=1:signal_length/fft_len

        y=W.*signal(n,:);                       % Output of channel estimator
        e = desired(n,:) - y;                   % Instantaneous error of LMS
        W = W + eta * e .* conj(signal(n,:));   % Weight update rule of LMS
        J(n) = e * e';                          % Instantaneuous squared error
        
    end
   temp1 = temp1 + W;
   temp2 = temp2 + J;
end
% Averaging results
temp1 = temp1./runs;  
temp2 = temp2./runs;
%% Results
% Impulse response
[bb]=invfreqz(conj(temp1'),Freq_vector,'complex',3,0);     % Estimated channel impulse response
[channel_impulse_response;bb]                    % Actual channel impulse response
% Estimation error and normalized impulse response difference
[norm(com_freq-conj(temp1')) 10*log10(norm(bb-channel_impulse_response)./norm(channel_impulse_response))] 

% Input and output signal estimation
r=ifft(sys_out,[],1);               % channel output in time domain
rnoisy=ifft(desired,[],1);          % noisy signal in time domain

model_out=signal.*repmat(temp1,signal_length/fft_len,1); % estimated channel output in frequency domain
mr=ifft(model_out,[],1);                             % estimated channel output in time domain

est_in=desired./repmat(temp1,signal_length/fft_len,1);   % estimated channel input in frequency domain
est_ir=ifft(est_in,[],1);                            % estimated channel input in time domain

% Mean squared error (MSE)
MSE = temp2./signal_length;        % Mean squared error of LMS

%% Plots
% Cost function plot
figure
fsize=14; % plot text font size
plot(10*log10(MSE),'','linewidth',4)
grid minor
xlabel('Iterations','FontName','Times New Roman','FontSize',fsize);
ylabel('Mean squared error (MSE) in (dB)','FontName','Times New Roman','FontSize',fsize);
title('Cost function (MSE vs epochs iteration)','FontName','Times New Roman','FontSize',6*fsize/5);
set(gca,'FontName','Times New Roman','FontSize',fsize)

% Signal Estimation
figure
% Estimated Input and Actualy input signal
subplot(1,2,1)
hold on
plot(reshape((est_ir),1,signal_length),'+r')
plot(reshape((x),1,signal_length),'o')
legend('estimated input','Actual input')
xlabel('Real axis','FontName','Times New Roman','FontSize',fsize);
ylabel('Imaginary axis','FontName','Times New Roman','FontSize',fsize);
title('constellation diagram','FontName','Times New Roman','FontSize',6*fsize/5);
set(gca,'FontName','Times New Roman','FontSize',fsize)

% Estimated Output and Actualy Output signal
subplot(1,2,2)
hold on
plot(reshape((r),1,signal_length),'+k')
plot(reshape((rnoisy),1,signal_length),'ok')
plot(reshape((mr),1,signal_length),'^r')
legend('Actual output','noisy output','estimated output')
xlabel('Real axis','FontName','Times New Roman','FontSize',fsize);
ylabel('Imaginary axis','FontName','Times New Roman','FontSize',fsize);
title('constellation diagram','FontName','Times New Roman','FontSize',6*fsize/5);
set(gca,'FontName','Times New Roman','FontSize',fsize)