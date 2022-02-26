function s_s = spectral_subtraction(y_fft, noise_est)

L = 5;
[~, fr_count] = size(y_fft);

s_s = zeros(size(y_fft));
s_s(:,1:L) = sqrt(max(1 - (noise_est(:,1:L)./y_fft(:,1:L)),0.2)) .* abs(y_fft(:,1:L));
for i = L : fr_count
    y_k = (1/L).*sum(abs(y_fft(:,i-L+1:i)).^2,2);
    est = sqrt(max(1 - (noise_est(:,i)./y_k),0.2));
    s_s(:,i) = est.* abs(y_fft(:,i));
end
end