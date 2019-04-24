function z = lms_equalizer(y, info, delta)
estimated_c = [0 0 0 0 0 1 0 0 0 0 0];
K = 5;
for k = 1 : size(y, 2) - 2*K
    y_k = y(k : k+2*K);
    
    z_k = estimated_c * y_k';
    e_k = info(k) - z_k;
    estimated_c = estimated_c +  delta * e_k * y_k;
    
    z(k) = z_k;
end
end