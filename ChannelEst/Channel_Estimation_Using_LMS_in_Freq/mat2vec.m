function result = mat2vec(mat)
% 将一个矩阵转换为行向量
% 接收参数为mat，期望操作是将mat按行拼接为一个行向量
% 函数输出一个转换后的行向量
    mat = mat.';
    result= mat(:).';
end