% 实现 Alamouti STBC 方案
function result = genSTBC(data)
stbc_temp = kron(reshape(data, 2, []), ones(1, 2));
stbc_temp(1 , 2 : 2 : end) = conj(stbc_temp(1 , 2 : 2 : end));
stbc_temp(2 , 2 : 2 : end) = -conj(stbc_temp(2 , 2 : 2 : end));
stbc_temp(:,   2 : 2 : end) = flip(stbc_temp(:, 2: 2: end));
result = stbc_temp;
end
