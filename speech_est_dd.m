function [speech_est, snr] = speech_est_dd(y_fft, noise_est)

[~, fr_count] = size(y_fft);

alpha_dd = 0.95;
speech_est = zeros(size(y_fft));
snr = zeros(size(y_fft));

speech_est(:,1) = abs(y_fft(:,1).^2) - noise_est(:,1);

for i = 1 : fr_count
    snr(:,i) = alpha_dd * (speech_est(:,i)./noise_est(:,i)) + (1 - alpha_dd)* max((abs(y_fft(:,i)).^2)./noise_est(:,i)-1,0);    
    speech_est(:,i) = snr(:,i).* noise_est(:,i);
end

end