function p_corr = sidak_corr(p, N)

p_corr = 1 - (1-p).^N;

end