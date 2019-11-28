function H_Est = refactor_LMS_Estimator(x, y)
%% Initializing simulation parameters
tempVar = 0;     % temporary estimated coefficients of channel (for each run)
pathNum = length(channelState);
       
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
%% LMS parameter
W = randn(size(com_freq'));                % Initial weights of LMS

    for n=1:signal_length/fft_len

        y=W.*signal(n, :);                       % Output of channel estimator
        e = desired(n, :) - y;                   % Instantaneous error of LMS
        W = W + eta * e .* conj(signal(n, :));   % Weight update rule of LMS
        J(n) = e * e';                          % Instantaneuous squared error       
    end
%% Results
% Impulse response
H_Est=invfreqz(conj(W'), Freq_vector, 'complex', pathNum - 1, 0);     % Estimated channel impulse response

end