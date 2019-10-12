function H_Est = LMSChannelEstimator(channelState)
%% Initializing simulation parameters
runs = 2;     % Number of Monte Carlo simulations (runs)
temp1 = 0;     % temporary estimated coefficients of channel (for each run)
N = (2^10);    % Number of samples
Bits = 2;      % For modulation    
fft_len=64;    % Fast Fourier Transform (Frame size) or channel length
eta = 1e-2;    % Learning rate for LMS
pathNum = length(channelState);

%% Defining Unknown channel
% Channel impulse response
channel_impulse_response = channelState;    
% Frequency response of channel
[com_freq, Freq_vector]=freqz(channel_impulse_response, 1, fft_len);

%% Monte Carlos Simulation
for k = 1 : runs
    % Generate Random signal for each independent run
    data = randi([0 ((2^Bits)-1)], 1, N);  % Random signal
    x = pskmod(data, 2^Bits);          % Phase shit keying (PSK) modulation
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

%% LMS parameter
W = randn(size(com_freq'));                % Initial weights of LMS

    for n=1:signal_length/fft_len

        y=W.*signal(n, :);                       % Output of channel estimator
        e = desired(n, :) - y;                   % Instantaneous error of LMS
        W = W + eta * e .* conj(signal(n, :));   % Weight update rule of LMS
        J(n) = e * e';                          % Instantaneuous squared error
        
    end
   temp1 = temp1 + W;
end
% Averaging results
temp1 = temp1./runs;  
%% Results
% Impulse response
H_Est=invfreqz(conj(temp1'), Freq_vector, 'complex', pathNum - 1, 0);     % Estimated channel impulse response


end