%% sample audio files
% load handel
% y = transpose(y);

% % sample audio file 2
% [y,Fs] = audioread('/MATLAB Drive/mixkit-melodical-flute-music-notification-2310.wav');

% % sample audio file 3
%  [y,Fs] = audioread('/MATLAB Drive/mix_07s_audio-joiner.com.wav');

y(35000:36000) = 0;
Y = y(30001:40000);

% sound(Y);
%% additive noise model
L = length(Y);  %length of speech signal 
speech = Y;
noise = 0.25*randn(1,length(speech)); % generating noise 
W = hann(L); % hanning window
W = transpose(W);
sk = speech.*W; % windowed speech
nk = noise.*W; % windowed noise 
xk =  sk + nk; % noise corrupted signal


% Xk = fft(xk);
[Xk,f,t] = stft(xk);
phase = angle(Xk);
% Nk = fft(nk);
Nk = stft(nk);
U = sum(abs(Nk))/length(nk);
M = 100;
X = abs(Xk);
N = length(X)-M+1;
n = 0:1/Fs:(length(Y)-1)/Fs;
%% plotting the additive noise model
figure(1);
sgtitle('The Additive Noise Model');
subplot(3,1,1)
plot(n,sk);
grid on;
xlabel('n');
title('Windowed Speech Signal');

subplot(3,1,2)
plot(n,nk,'blue');
grid on;
xlabel('n');
title('Windowed Noise Signal');

subplot(3,1,3)
plot(n,xk, 'm');
grid on;
grid on;
xlabel('n');
title('Windowed Noise Corrupted Signal');

%% Noisy signal Spectrogram : 
figure(3)
imagesc(t,f,(abs(Xk))),axis xy,xlabel('t(secondes)'),ylabel('frequency(Hz)')
title(' Spectrogram of "Mix" ( signal + noise )'),colormap jet

%% Magnitude averaging
avgX = zeros(1,N);

% taking the local average 
for k = 1:length(Xk)-M+1
    avgX(k) = sum(X(k:k+M-1))/M;
    
end

len = length(avgX);
w = 0:(2*pi)/len:(2*pi*(len-1))/len;
phase1 = phase(1:length(avgX));
S = (abs(avgX)-U(1:length(avgX))).*exp(1j*phase1);

%% plotting fft of the windowed speech signal 
figure(2);
subplot(2,2,1);
imagesc(t,f,(abs(S))),axis xy,xlabel('t(secondes)'),ylabel('frequency(Hz)')
title('Spectrogram after magnitude averaging'),colormap jet

%% Half wave rectification
H = 1 - U(1:length(avgX))./avgX; 
Hr = (H + abs(H))/2;
S1 = Hr.*avgX; % estimate after half wave rectification

%% plotting spectral subtraction estimator after HWR
subplot(2,2,2);
imagesc(t,f,(abs(S1))),axis xy,xlabel('t(secondes)'),ylabel('frequency(Hz)')
title('Spectrogram after HWR'),colormap jet

%% Residual Noise Reduction
NR = Nk - U.*exp(angle(Nk)); % residual noise
max_NR = max(NR); % max residual noise
S2 = zeros(size(S1));

for k = 1:length(avgX)
    
    % taking the minimum estimate from 3 adjacent frames 
    if(abs(S1(k))<max_NR) 
        Nframe = [abs(S1(k-1)),abs(S1(k)),abs(S1(k+1))];
        S2(k) = min(Nframe);
    else
        S2(k) = S1(k);
    end
end

%% plotting spectral subtraction estimator after RNR
subplot(2,2,3);
imagesc(t,f,(abs(S2))),axis xy,xlabel('t(secondes)'),ylabel('frequency(Hz)')
title('Spectrogram after RNR'),colormap jet

%% Additional Signal Attenuation
I = abs(S2)./abs(U(1:length(avgX))); 
T = 20*log10(I); % indicates the absense of speech (SNR value)
S3 = zeros(size(S2));

c = 10^(-3/2);
for k = 1:length(avgX)
    % absense of speech
    if(T(k) <= -20)  
        S3(k) = c*Xk(k);
    else 
        S3(k) = S2(k);
    end
end
%% plotting the speech estimater after ASA
subplot(2,2,4);
imagesc(t,f,log(abs(S3))),axis xy,xlabel('t(secondes)'),ylabel('frequency(Hz)')
title('Spectrogram after magnitude averaging'),colormap jet

%% plotting after each error reduction in frequency domain
% plotting spectral subtraction estimator after Magnitude averaging
figure(6);
subplot(2,2,1);
plot(w,abs(S(1:length(avgX))));
grid on;
title('Speech spectrum after magnitude averaging');
xlabel('w');
ylabel('S[e^jw]');

% plotting spectral subtraction estimator after Half wave rectification
subplot(2,2,2);
plot(w,abs(S1(1:length(avgX))));
grid on;
title('Speech spectrum after Half wave rectification');
xlabel('w');
ylabel('S[e^jw]');

% plotting spectral subtraction estimator after Residual noise removal
subplot(2,2,3);
plot(w,abs(S2(1:length(avgX))));
grid on;
title('Speech spectrum after residual noise removal');
xlabel('w');
ylabel('S[e^jw]');

% plotting after additional signal attenuation
subplot(2,2,4);
plot(w,abs(S3(1:len)));
grid on;
title('Speech spectrum after additional signal attenuation');
xlabel('w');
ylabel('S[e^jw]');
%% reconstructing the noise reduced signal 
Output = ifft(S3);
figure(5);
imagesc(t,f,abs(Output)),axis xy,xlabel('t(secondes)'),ylabel('frequency(Hz)')
title('Spectrogram after spectral subtraction'),colormap jet
% sound(real(Output),Fs);
