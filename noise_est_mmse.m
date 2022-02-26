function psd_n3 = noise_est_mmse(y_psd, y_fft)

[~ , fr_count] = size(y_psd);

psd_n3 = zeros(size(y_psd));
exp_noise = zeros(size(y_psd));
alpha_mmse = 0.85;
snr_h1 = 10^1.5;

% Prior Probabilities of Speech and Noise Presence
p_0 = 0.5;
p_1 = 0.5;

psd_n3(:,1) = abs(y_fft(:,1).^2);

for i = 2 : fr_count
        p_y_0 = (1./(pi*psd_n3(:,i-1))).*exp(-abs(y_fft(:,i).^2)./psd_n3(:,i-1)); 
        p_y_0(p_y_0 < 0.001) = 0.001;

        p_y_1 = (1./(pi*psd_n3(:,i-1).*(1+snr_h1))).*exp(-abs(y_fft(:,i).^2)./(psd_n3(:,i-1).*(1+snr_h1))); 
        p_y_1(p_y_1 < 0.001) = 0.001;      
        
%         temp = (1 + (p_1/p_0)*(1 + snr_h1)*exp((-abs(y_fft(:,i).^2)*snr_h1)./(psd_n3(:,i-1)*(1+snr_h1))));
%         p_h1 = 1/temp;
        p_h1 = (p_1 .*p_y_1)./((p_1.* p_y_1) + (p_0.* p_y_0));
        p_h0 = 1 - p_h1;
        
    exp_noise(:,i) = p_h0 .* abs(y_fft(:,i)).^2 + p_h1.*psd_n3(:,i-1);
    
    psd_n3(:,i) = alpha_mmse * psd_n3(:,i-1) + (1 - alpha_mmse)*(exp_noise(:,i));
end
end