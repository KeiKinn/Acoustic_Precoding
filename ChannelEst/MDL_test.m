clc
clear all
close all

% N = 10;
% d = 0.5;
% elementPos = (0:N-1)*d;
% angles = [0 -25];
% x = sensorsig(elementPos,300,angles,db2pow(-5));
% Nsig = mdltest(x)

baseData = importdata('Rk-noisy.mat');
data = zeros(7, length(baseData));
for counter = 1 : 7
    data(counter, counter : end) = baseData(1 : end - counter + 1);
end

anti_diag = flip(eye(7, 7));

R_FB = 0.5 * (data * data' + anti_diag * conj(data * data') * anti_diag);

eigb = sort(eig(R_FB), 'descend');

nsig = mdltest(eigb);