function g_cal = mmse_gain(speech, noise, y_ffta)
% [fr_size, fr_count] = size(speech);

p_beta = 1;
p_gamma = 1;
p_nu = 0.6;

% g_cal = zeros(size(speech));


% for k = 1 : fr_count
%     for j = 1 : fr_size
        
        p_A = @(a)(((p_gamma * p_beta.^p_nu)/(gamma(p_nu)))*(a.^(p_gamma * (p_nu - 1))) .*exp(-p_beta *a.^p_gamma));
        p_a_phi = @(a, phi)(a.*exp(1i*phi));
%        p_y_s = @(a, phi)((1/(pi*noise(j,k)))*exp((-abs(y_fft(j,k) - (a.*exp(1i*phi))))./noise(j,k)));
        p_y_s = @(a, phi)((1/(pi*noise))*exp((-abs(y_ffta - (a.*exp(1i*phi))))./noise));
        
        numer = @(a, phi) p_A(a) .* p_y_s(a, phi) .* p_a_phi(a,phi);
        denom = @(a, phi) p_A(a) .* p_y_s(a, phi);

%        g_cal(j,k) = integral2(numer, 0, inf, 0, 2*pi)/integral2(denom, 0, inf, 0, 2*pi);
        g_cal = integral2(numer, 0, inf, 0, 2*pi)/integral2(denom, 0, inf, 0, 2*pi);
%    end
end
