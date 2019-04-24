function [info] = random_binary(N)
if nargin == 0
    N = 1e4;
end

for i = 1 : N
    temp = rand;
    if (temp < 0.5)
        info(i) = -1;
    else
        info(i) = 1;
    end
end
end