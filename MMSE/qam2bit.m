% 将接收到的QAM转化为bit，以便于BER
function result = qam2bit(input, QAM_order)
rx_zf_qam_demod(1, :) = qamdemod(input(1, :), QAM_order);
rx_zf_qam_demod(2, :) = qamdemod(input(2, :), QAM_order);
rx_decode = zeros(1, 2 * length(rx_zf_qam_demod));
rx_decode(1 : 2 :end) = rx_zf_qam_demod(1, :);
rx_decode(2 : 2 :end) = rx_zf_qam_demod(2, :);
result = de2bi(rx_decode);
end