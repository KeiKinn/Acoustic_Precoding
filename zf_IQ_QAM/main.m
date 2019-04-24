clear;
clc;
echo off;
close all;

N = 10000;
info = random_binary(N);
SNR_in_dB = 8 : 1: 18;

for j = 1 : length(SNR_in_dB)
     zf_info = zf_precoding(info);
    [y, len] = channel(zf_info, SNR_in_dB(j));   
    numoferr = 0;
    for i = 1: N 
        if(y(i) < 0)
            decis = -1;
        else
            decis = 1;
        end
        
        if (decis ~= info(i))
            numoferr = numoferr + 1;
        end
    end
    
    Pe(j) = numoferr / N;
end

semilogy(SNR_in_dB, Pe, 'red*-');
hold on;

delta_1 = 0.11;
delta_2 = 0.09;

for j = 1 : length(SNR_in_dB)
    y = channel(info, SNR_in_dB(j));
    z= lms_equalizer(y, info, delta_1);
    numoferr = 0;
    for i = 1 : N
        if(z(i) <0)
            decis = -1;
        else
            decis = 1;
        end
        
        if(decis ~= info(i))
            numoferr = numoferr + 1;
        end
    end
    Pe(j) = numoferr / N;
end

semilogy(SNR_in_dB, Pe, 'blacko-');
hold on;

for j = 1 : length(SNR_in_dB)
    y = channel(info, SNR_in_dB(j));
    z= lms_equalizer(y, info, delta_2);
    numoferr = 0;
    for i = 1 : N
        if(z(i) <0)
            decis = -1;
        else
            decis = 1;
        end
        
        if(decis ~= info(i))
            numoferr = numoferr + 1;
        end
    end
    Pe(j) = numoferr / N;
end

semilogy(SNR_in_dB, Pe, 'blue.-');
hold off;
figloc

legend('no equalizer', 'delta = 0.11', 'delta = 0.09 ');