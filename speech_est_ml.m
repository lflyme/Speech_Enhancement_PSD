function [speech_est,snr] = speech_est_ml(y_fft, noise_est)

K = 20;
[~, fr_count] = size(y_fft);

speech_est = zeros(size(y_fft));
snr = zeros(size(y_fft));

speech_est(:,1) = abs(y_fft(:,1).^2) - noise_est(:,1);
snr(:,1:K) = max((abs(y_fft(:,1:K)).^2./noise_est(:,1:K) )  - 1,0);


for i = K : fr_count
    numer = (1/K).*sum(abs(y_fft(:,i-K+1:i)).^2,2);
    denom = noise_est(:,i);
    snr(:,i) = max((numer ./ denom - 1),0);
    speech_est(:,i) = snr(:,i).* noise_est(:,i);
end

end