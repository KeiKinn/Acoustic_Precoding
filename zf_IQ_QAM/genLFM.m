function LFMData = genLFM(f_start, f_end, time, fs)
f1 = f_start;
f2 = f_end;
t = 0 : 1 / fs : (time - 1 / fs);
M = (f2 - f1) / max(t);
LFMData = cos(2 * pi * (f1 * t + 0.5 * M * t.^2));
end