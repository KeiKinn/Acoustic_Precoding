% 为2x2的STBC系统求解信道补偿矩阵
% 在本系统中，因为信号首先经过了预编码矩阵W
% 因此本身对H求伪逆就变成了对等效信道矩阵HW求伪逆
function result = invPseH(pseChannel)
H_base = pseChannel;
H_fliplr = conj(fliplr(H_base));
H_fliplr(:, 2) = -H_fliplr(:, 2);
H = [H_base; H_fliplr];
result = (H' * H) \ H';
end