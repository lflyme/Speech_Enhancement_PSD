function [psd_n1] = noise_est_ms(y_psd,y_fft)

[fr_size, fr_count] = size(y_psd);

D = 20;
M = 0.705;

fs = 16000;
fr_overlap = 256;
L = fr_size;
R = fr_overlap;
T_SM = 0.2;
alpha_ms = (T_SM *fs/R - 1)/ (T_SM *fs/R + 1);
beta_ms = 0.85;

alpha_min = min(0.3, 15^(-R/(0.064 * fs)));

psd_n1 = zeros(size(y_psd));
p_yy = zeros(size(y_psd));
p_min = zeros(size(y_psd));

p_yy(:,1) = (2*pi*L).*abs(y_fft(:,1)).^2;
psd_n1(:,1) = y_psd(:,1);

p_fm = zeros(size(y_psd));
p_sm = zeros(size(y_psd));

p_fm(:,1) = (2*pi*L).*abs(y_fft(:,1)).^2;
p_sm(:,1) =(2*pi*L).*abs(y_fft(:,1)).^4;

for k = 2 : fr_count  
	% A posteriori SNR
    snr_posteriori =  ((2*pi*L).*abs(y_fft(:,k-1)).^2) ./ psd_n1(:,k-1);
    
    % PSD Smoothening Parameter for PSD Update
    alpha_ms = 1./((1 + (snr_posteriori - 1).^2));
    if(alpha_ms > 0.96.*ones(fr_size,1))
        alpha_ms = 0.96.*ones(fr_size,1);
    end 
%     if(alpha_ms < alpha_min.*ones(fr_size,1))
%         alpha_ms = alpha_min.*ones(fr_size,1);
%     end
    
	% Smoothed Periodogram
    p_yy(:,k) = alpha_ms .*(p_yy(:,k-1)) + (1 - alpha_ms).*abs(y_fft(:,k).^2);  
    
    if(k > D)
        p_min(:,k) = min(p_yy(:,k-D+1:k),[],2);
    else
        p_min(:,k) = min(p_yy(:,1:k),[],2);
    end
   
    % PSD Smoothening Parameter for Moments
    beta_ms = alpha_ms.^2;
    if(beta_ms > 0.8.*ones(fr_size,1))
        beta_ms = 0.8.*ones(fr_size,1);
    end
        
    % First and Second Moments of PSD Estimate
    p_fm(:,k)  = (beta_ms.* p_fm(:,k-1)) + ((1 - beta_ms).* y_psd(:,k));
    p_sm(:,k)  = (beta_ms.* p_sm(:,k-1)) + ((1 - beta_ms).* (y_psd(:,k).^2));
    
    % Variance of PSD Estimate
    var_p_s = p_sm(:,k) - (p_fm(:,k).^2);
    
    % Equivalent Degree of Freedom
    Q_eq = (2 .* psd_n1(:,k-1).^2)./ var_p_s ;
    
    % Equivalent Degree of Freedom - Scaled
    Q_eq_h = (Q_eq - (2 * M) )./(1 - M);
    
    % Bias Compensation
    B_min = 1 + ((2.*(D - 1))./Q_eq_h);
    
    Q_inv = (1/L) .* sum(1./Q_eq);
    B_c = 1 + 2.12 *sqrt(Q_inv);
    
    psd_n1(:,k) = p_min(:,k) .* B_min; % .* B_c;        
end

end