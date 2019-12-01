function [result, randi_bit] = genQAM(order, symbolPerTime)
% 生成 QAM，如果不指定参数，则仅默认生成16QAM，1000个码元
if nargin == 1
    symbolPerTime = 1000;
end
if nargin == 0
    order = 16;
    symbolPerTime = 1000;
end

k = log2(order);
randi_dec = randi(order, symbolPerTime, 1) - 1;
randi_bit = de2bi(randi_dec, k);
result = qammod(randi_dec, order);
end