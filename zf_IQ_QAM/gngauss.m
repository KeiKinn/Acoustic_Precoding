function [gsrv1, gsrv2] = gngauss(m, sgma)
if nargin == 0
    m = 0;
    sgma = 1;
elseif nargin == 1
    sgma = m;
    m = 0;
end

u =rand;
z = sgma * (sqrt(2 * log(1 / (1 - u))));
u =rand;
gsrv1 = m + z * cos(2 *pi * u);
gsrv2 = m + z * sin(2 * pi * u);
end