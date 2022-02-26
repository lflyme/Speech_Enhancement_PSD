function [s_w] = weiner_gain(snr_est, y_fft)
g_weiner = snr_est./(1 + snr_est);
g_weiner(g_weiner < 0.1) = 0.1; 
s_w = g_weiner .* y_fft;
end