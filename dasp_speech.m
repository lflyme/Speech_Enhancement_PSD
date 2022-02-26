warning off
clc;
clear all;
close all;

[speech, fs]   = audioread('clean_speech.wav');
[noise_1, ~] = audioread('aritificial_nonstat_noise.wav');
[noise_2, ~] = audioread('babble_noise.wav');
[noise_3, ~] = audioread('Speech_shaped_noise.wav');

SNR = 15;

% modified for clarity
speech_norm = sqrt(mean(speech.^2)/(10^(SNR/10)));

noise_1 = speech_norm *noise_1(1:length(speech))./mean(mean(noise_1(1:length(speech)).^2)).^(0.5);
noise_2 = speech_norm *noise_2(1:length(speech))./mean(mean(noise_2(1:length(speech)).^2)).^(0.5);
noise_3 = speech_norm *noise_3(1:length(speech))./mean(mean(noise_3(1:length(speech)).^2)).^(0.5);
noise_4 = speech_norm *randn(size(speech));

noise_n = noise_4;

y = speech + noise_n;
s = speech;
n = noise_n;

%% Segmentation : Signals to Matrix of Frames

t_seg = 32e-3;
fr_size  = t_seg * fs;
fr_overlap = 0.5 * fr_size;
fr_count = floor(length(y)/fr_overlap);

y = y(1:fr_count*fr_overlap);
s = s(1:fr_count*fr_overlap);
n = n(1:fr_count*fr_overlap);

% Frame Initializations
fr_unit = zeros(fr_size,1);
y_frames = zeros(fr_size,1);
s_frames = zeros(fr_size,1);
n_frames = zeros(fr_size,1);


for i = 1 : fr_overlap : ((fr_count -2) * fr_overlap) + 1
    % Frame selection - Current index to next (fr_size) indices
    fr_unit = i : i + fr_size - 1;
    % Signal Frames
    y_frames = cat(2,y_frames,y(fr_unit));
    s_frames = cat(2,s_frames,s(fr_unit));
    n_frames = cat(2,n_frames,n(fr_unit));
end

y_frames(:,1) = [];
s_frames(:,1) = [];
n_frames(:,1) = [];

%% Discrete Fourier Transform
% FFT of Signals
y_fft = fft(y_frames);
s_fft = fft(s_frames);
n_fft = fft(n_frames);

%% Power Spectral Density
% PSD for Signals
L = fr_size;
D = 80;
y_psd = 2*pi*L*periodogram(y_frames,rectwin(L),L,'twosided');
s_psd = 2*pi*L*periodogram(s_frames,rectwin(L),L,'twosided');
n_psd = 2*pi*L*periodogram(n_frames,rectwin(L),L,'twosided');

%% Noise Estimation Methods

% A. Minimum Statistics Estimation
[noise_ms] = noise_est_ms(y_psd,y_fft);

% B. Voice Activity Detector Based Estimation
[noise_vad] = noise_est_vad(y_psd, y_frames, y_fft);

% C. MMSE Based Estimation
[noise_mmse] = noise_est_mmse(y_psd, y_fft);

%% Speech Estimation Methods

% A. Maximum Likelihood 
[speech_est_ml_mmse, snr_ml_mmse] = speech_est_ml(y_fft, noise_mmse);
[speech_est_ml_ms, snr_ml_ms] = speech_est_ml(y_fft, noise_ms);
[speech_est_ml_vad, snr_ml_vad] = speech_est_ml(y_fft, noise_vad);

% B. Decision Directed Approach
[speech_est_dd_mmse, snr_dd_mmse] = speech_est_dd(y_fft, noise_mmse);
[speech_est_dd_ms, snr_dd_ms] = speech_est_dd(y_fft, noise_ms);
[speech_est_dd_vad, snr_dd_vad] = speech_est_dd(y_fft, noise_vad);


%% Target Gain Estimation Methods

% A. Weiner Smoother Gain Heuristic
s_w_dd_ms = weiner_gain(snr_dd_ms, y_fft);
s_w_dd_vad = weiner_gain(snr_dd_vad, y_fft);
s_w_dd_mmse = weiner_gain(snr_dd_mmse, y_fft);

s_w_ml_ms = weiner_gain(snr_ml_ms, y_fft);
s_w_ml_vad = weiner_gain(snr_ml_vad, y_fft);
s_w_ml_mmse = weiner_gain(snr_ml_mmse, y_fft);


% B. Spectral Subtraction Estimation
s_s_ms = spectral_subtraction(y_fft, noise_ms);
s_s_vad = spectral_subtraction(y_fft, noise_vad);
s_s_mmse = spectral_subtraction(y_fft, noise_mmse);


% % C. MMSE Heuristic
% g_m = zeros(size(y_psd));
% for k = 1 : fr_count
%     for j = 1 : fr_size
%         g_m(j,k) = mmse_gain(speech_dd(k,j), noise_mmse(k,j),y_fft(k,j));
%     end
% end     
% s_mmse = g_m .* y_fft;


%% Post-Processing of Frames

% Inverse Fourier Transform on Frames + Overlap and Addition for Signal
% Reconstruction

% Results from Decision Directed Approach
s_w_p_dd_ms = real(ifft(s_w_dd_ms));
speech_w_dd_ms = overlap_add(s_w_p_dd_ms);

s_w_p_dd_vad = real(ifft(s_w_dd_vad));
speech_w_dd_vad = overlap_add(s_w_p_dd_vad);

s_w_p_dd_mmse = real(ifft(s_w_dd_mmse));
speech_w_dd_mmse = overlap_add(s_w_p_dd_mmse);

% Results from Maximum Likelihood Approach
s_w_p_ml_ms = real(ifft(s_w_ml_ms));
speech_w_ml_ms = overlap_add(s_w_p_ml_ms);

s_w_p_ml_vad = real(ifft(s_w_ml_vad));
speech_w_ml_vad = overlap_add(s_w_p_ml_vad);

s_w_p_ml_mmse = real(ifft(s_w_ml_mmse));
speech_w_ml_mmse = overlap_add(s_w_p_ml_mmse);

% Results from 
s_s_p_ms = real(ifft(s_s_ms));
speech_s_ms= overlap_add(s_s_p_ms);

s_s_p_vad = real(ifft(s_s_vad));
speech_s_vad = overlap_add(s_s_p_vad);

s_s_p_mmse = real(ifft(s_s_mmse));
speech_s_mmse = overlap_add(s_s_p_mmse);


%% Evaluations of the Methods

speech = speech(1:length(speech_w_dd_mmse));
% Mean Squared Errors
mse_dd_mmse = mean((speech - speech_w_dd_mmse).^2);
mse_dd_ms = mean((speech - speech_w_dd_ms).^2);
mse_dd_vad = mean((speech - speech_w_dd_vad).^2);


mse_ml_mmse = mean((speech - speech_w_ml_mmse).^2);
mse_ml_ms = mean((speech - speech_w_ml_ms).^2);
mse_ml_vad = mean((speech - speech_w_ml_vad).^2);


mse_s_mmse  = mean((speech - speech_s_mmse).^2);
mse_s_ms  = mean((speech - speech_s_ms).^2);
mse_s_vad  = mean((speech - speech_s_vad).^2);

mse = [mse_dd_ms, mse_dd_vad, mse_dd_mmse; mse_ml_ms, mse_ml_vad, mse_ml_mmse; mse_s_ms,  mse_s_vad,  mse_s_mmse];

% Mean Absolute Errors
mae_dd_mmse = mean((speech - speech_w_dd_mmse));
mae_dd_ms = mean((speech - speech_w_dd_ms));
mae_dd_vad = mean((speech - speech_w_dd_vad));


mae_ml_mmse = mean((speech - speech_w_ml_mmse));
mae_ml_ms = mean((speech - speech_w_ml_ms));
mae_ml_vad = mean((speech - speech_w_ml_vad));


mae_s_mmse  = mean((speech - speech_s_mmse));
mae_s_ms  = mean((speech - speech_s_ms));
mae_s_vad  = mean((speech - speech_s_vad));

mae = [mae_dd_ms, mae_dd_vad, mae_dd_mmse; mae_ml_ms, mae_ml_vad, mae_ml_mmse; mae_s_ms, mae_s_vad, mae_s_mmse];


% SNR Comparison
noise_n = noise_n(1:length(speech));

snr_orig = snr(speech,noise_n);

snr_dd_ms = snr(speech_w_dd_ms,noise_n); 
snr_dd_vad = snr(speech_w_dd_vad,noise_n); 
snr_dd_mmse = snr(speech_w_dd_mmse,noise_n); 

snr_ml_ms = snr(speech_w_ml_ms,noise_n); 
snr_ml_vad = snr(speech_w_ml_vad,noise_n); 
snr_ml_mmse = snr(speech_w_ml_mmse,noise_n); 

snr_s_ms = snr(speech_s_ms,noise_n); 
snr_s_vad = snr(speech_s_vad,noise_n); 
snr_s_mmse = snr(speech_s_mmse,noise_n); 

snrs = [snr_dd_ms, snr_dd_vad, snr_dd_mmse; snr_ml_ms, snr_ml_vad, snr_ml_mmse; snr_s_ms, snr_s_vad, snr_s_mmse];

%% Write Output Audio

audiowrite('noisy.wav' , y/max(abs(y)), fs);
audiowrite('enhanced_w_dd_mmse.wav' , speech_w_dd_mmse/max(abs(speech_w_dd_mmse)), fs);

%% Plots

a = figure(1);
subplot(3,1,1);
plot(speech)
title('Clean speech')
subplot(3,1,2);
plot(y)
title('Noisy speech')
subplot(3,1,3);
plot(speech_w_dd_mmse)
title('Enhanced speech')
saveas(a,'over.png');


b = figure(2);
plot(10*log(noise_ms(104,1:500)))
hold on;
plot(10*log(noise_vad(104,1:500)))
hold on;
plot(10*log(noise_mmse(104,1:500)))
hold on;
plot(10*log(n_psd(104,1:500)))
title('Estimated Noise \sigma_n^2 with Different Approaches')
legend("\sigma_nMS^2","\sigma_nVAD^2","\sigma_nMMSE^2","Periodogram", 'Location', 'best');
saveas(b,'noise_comp.png');


c = figure(3);
plot(10*log(speech_est_dd_ms(104,1:500)))
hold on;
plot(10*log(speech_est_dd_vad(104,1:500)))
hold on;
plot(10*log(speech_est_dd_mmse(104,1:500)))
hold on;
plot(10*log(s_psd(104,1:500)))
title('Estimated Noise \sigma_s^2 with Different Approaches with DD Approach')
legend("\sigma_sMS^2","\sigma_sVAD^2","\sigma_sMMSE^2","Periodogram",'Location', 'best');
saveas(c,'speech_dd.png');

d = figure(4);
plot(10*log(speech_est_ml_ms(104,1:500)))
hold on;
plot(10*log(speech_est_ml_vad(104,1:500)))
hold on;
plot(10*log(speech_est_ml_mmse(104,1:500)))
hold on;
plot(10*log(s_psd(104,1:500)))
title('Estimated Noise \sigma_s^2 with Different Approaches with ML Approach')
legend("\sigma_sMS^2","\sigma_sVAD^2","\sigma_sMMSE^2","Periodogram", 'Location', 'best');
saveas(d,'noise_ml.png');

approach = {"Weiner + DD","Weiner + DD","Weiner + DD"; "Weiner + ML","Weiner + ML","Weiner + ML"; "Spectral Subtraction", "Spectral Subtraction" ,"Spectral Subtraction"  };
e = figure(5);
bar(mse)
set(gca,'xticklabel',approach);
title('MSE for Various Noise Estimation Approahces')
xlabel('Methods')
ylabel('MSE')
legend('MS', 'VAD', 'MMSE','Location', 'best');
saveas(e,'mse.png');

f = figure(6);
bar(mae)
set(gca,'xticklabel',approach);
title('MAE for Various Noise Estimation Approahces')
xlabel('Methods')
ylabel('MAE')
legend('MS', 'VAD', 'MMSE','Location', 'best');
saveas(f,'mae.png');

g = figure(7);
bar(snrs)
yline(snr_orig,'-','True SNR')
set(gca,'xticklabel',approach);
title('SNR for Various Noise Estimation Approahces')
xlabel('Methods')
ylabel('SNR')
legend('MS', 'VAD', 'MMSE','Location', 'best');
saveas(g,'snr.png');