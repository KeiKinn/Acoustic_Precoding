function [channel_gain, channel_delay] = genChannelop(channel_order)
% 从张歆老师程序DSSS_FSK5.m中找到的3个信道，编码为1，2，3
channel_order = 1;
switch channel_order
    case 1
        rpt =[134.227    65.569
                  134.879   66.757
                  137.526   66.396
                  142.250   65.923];
    case 2
        rpt =[134.2152   65.867
                 140.5465   66.061
                 135.7882   65.906
                 134.5818   65.971
                 135.0426   65.951
                 136.2212   66.079
                 141.2072   66.315];
    case 3
        rpt = [121.2657   46.1790
                  121.1627   46.3129
                  126.2940   46.1332
                  128.4948   46.2669
                  143.3866   46.1812
                  128.2482   46.4922
                  132.3180   46.4071
                  128.5326   46.2748
                  129.2960   46.2072
                  128.8532   46.3496
                  131.8432   46.3626
                  128.9838   46.3262
                  133.0679   46.3668
                  141.5868   46.2857];
end
channel_gain = rpt(:,1);
channel_delay = rpt(:,2);
channel_gain = min(channel_gain) - channel_gain;
for j =1:length(channel_gain)
    channel_gain(j) = 10^(channel_gain(j) / 20);
end
channel_delay = channel_delay - min(channel_delay);
        
end