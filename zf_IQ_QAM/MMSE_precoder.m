function result = MMSEPrecoder(H)
[U, ~] = size(H);
result = H' / (H*H');
end