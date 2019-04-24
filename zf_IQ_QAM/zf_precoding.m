function coded_data = zf_precoding(data)
actual_isi = [0.05 -0.063 0.088 -0.126 -0.25 0.9047 0.25 0 0.126 0.038 0.088];
w_len =length(-7:7);
P = convmtx(actual_isi', w_len);
u_ZF = zeros(25, 1);
u_ZF(13) = 1;
c_Ls = ((P'*P))\((P'))*u_ZF;

coded_temp = filter(c_Ls, 1, data);
coded_data = coded_temp(11 : end);
end