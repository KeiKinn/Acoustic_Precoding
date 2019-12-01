% 将 2x2 转换为 4x1 以便于STBC解码
% input：2 x N（N>2）
function result = shapeMat(input)
result= input(:, 1 : 2 : end);
result = [result; conj(input(:, 2 : 2 : end))];
end