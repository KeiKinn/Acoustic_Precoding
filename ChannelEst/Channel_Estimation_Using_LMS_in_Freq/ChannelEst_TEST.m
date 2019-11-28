clear all;
close all;
clc


%% - - - Channel- - - %%
H_base = [0.8, 0.5 + 0.1j; 0.3 + 0.2j, 0.6];
[~, col] = size(H_base);
H_Est = zeros(col);
for counter_i = 1 : col
    H_Est(counter_i, :) = LMSChannelEstimator(H_base(counter_i, :));
end