function H_Est = LSEstimator(channelState, signalForm, MentoCarloNum)
if nargin < 3
    MentoCarloNum = 10; % Number of Monte Carlo simulations (runs)
end

%% Initializing simulation parameters
tempVar = 0;     % temporary estimated coefficients of channel (for each run)
pathNum = length(channelState);

%% QPSK Para
N = (2^10);    % Number of samples
Bits = 2;      % For modulation
fft_len=64;    % Fast Fourier Transform (Frame size) or channel length
eta = 0.2e-2;    % Learning rate for LMS

%% Defining Unknown channel
channel_impulse_response = channelState;    % Channel impulse response
[com_freq, Freq_vector]=freqz(channel_impulse_response, 1, fft_len); % Frequency response of channel

%% Monte Carlos Simulation

H_Est_Temp = [];
for k = 1 : MentoCarloNum
    if signalForm == 0
        % Generate Random signal for each independent run
        data = randi([0 ((2^Bits)-1)], 1, N);  % Random signal
        x = pskmod(data, 2^Bits);          % Phase shit keying (PSK) modulation
    else
        x =  lteZadoffChuSeq(1, N - 1)'; % ZC sequnence
    end
       
    % zero padding to ensure each frame is of size equals to the channel/frame
    padding=rem(length(x), fft_len);
    if padding~=0
        padding=fft_len-padding;
    end
    x=[x zeros(1, padding)];
    signal_length=length(x); % final length after padding
    
    % Dividing data stream into frames of length equals to fft_len
    x=reshape(x, [signal_length/fft_len fft_len]);
    
    % Convert time domain signal to frequency domain signal using FFT
    signal=fft(x, [], 1);
    
    % Generate desired output signal by (Convolution in time domain or
    % multiplication in frequency domain)
    sys_out=signal.*repmat(conj(com_freq'), signal_length/fft_len, 1);    % signal after passing through channel
    
    desired = sys_out;    % Addition of white gaussian noise (desired signal after disturbance)
%     desired = awgn(sys_out, SNR);
    
    %% LS parameter
for counter_i = 1 : signal_length / fft_len
   H_fre_temp = desired(counter_i, :)./signal(counter_i, :);
   H_Est_Temp = H_Est_Temp + H_fre_temp;
end
   tempVar = tempVar + H_Est_Temp / counter_i;
end
% Averaging results
tempVar = tempVar./MentoCarloNum;
%% Results
% Impulse response
H_Est=invfreqz(conj(tempVar'), Freq_vector, 'complex', pathNum - 1, 0);     % Estimated channel impulse response

end