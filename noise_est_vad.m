function psd_n2 = noise_est_vad(y_psd, y_frames, y_fft)

[fr_size , fr_count] = size(y_psd);

alpha_vad = 0.95;
kappa_vad = 0.9;

psd_n2 = zeros(size(y_psd));
snr_update = zeros(size(y_psd));
X_lap = zeros(size(y_psd));
gamma_lr_lap =  zeros(size(y_psd));

X_r = real(y_fft);
X_i = imag(y_fft);

psd_n2(:,1) = abs(y_fft(:,1).^2);
X_lap(:,1)  = abs(y_fft(:,1).^2);

for i = 2 : fr_count
    
    snr_aposteriori = (abs(y_fft(:,i).^2)./psd_n2(:,i-1));
    snr_apriori     = (abs(X_lap(:,i-1)).^2 ./psd_n2(:,i-1));
       
    X_lap(:,i) = abs( kappa_vad * abs(y_fft(:,i)) + (1 - kappa_vad) *sqrt(psd_n2(:,i-1) + abs(y_fft(:,i).^2)));
    
    snr_update(:,i) = alpha_vad * (snr_apriori) + (1 - alpha_vad)*(snr_aposteriori - 1);
    
    A = abs(X_r(:,i)) + abs(X_i(:,i));
    B = (X_lap(:,i) - sqrt(psd_n2(:,i-1)))./(X_lap(:,i) .* sqrt(psd_n2(:,i-1)));
    
    gamma_lr_lap(:,i) = (1./ sqrt(1 + snr_update(:,i))).* exp(sqrt(2/3).* A.* B); 
    
    
    for j = 1 : fr_size
        if(gamma_lr_lap(j,i) < 1)
            psd_n2(j,i) = alpha_vad .* psd_n2(j,i-1) + (1 - alpha_vad).* abs(y_fft(j,i-1).^2);
        else
            psd_n2(j,i) =  psd_n2(j,i-1);
        end
    end
end
