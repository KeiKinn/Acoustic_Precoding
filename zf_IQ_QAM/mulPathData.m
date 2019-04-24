function rx_data = mulPath(data, time_delay, channel_gain)
rx_data = zeros(1, length(data) + max(time_delay(1, :)));
for counter_i = 1 : length(time_delay)
    pre_zeros   = zeros(1, time_delay(1, counter_i));
    post_zeros = zeros(1, time_delay(2, counter_i));
    signal_temp = [pre_zeros, data, post_zeros];
    signal_temp = channel_gain(counter_i) * signal_temp;
    rx_data = rx_data + signal_temp;
end
end