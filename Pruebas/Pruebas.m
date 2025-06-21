%A = [1 -1/2 0 0 0 0 0 0 0 0 0 0]; 
% B = CalcularInversaZ(A); 
% subplot(3,1,1) 
% stem(A) 
% subplot(3,1,2) 
% stem(B) 
% subplot(3,1,3) 
% stem(conv(A,B))
clc
clear
Fs = 100; %100hz 
N = 5000; %Genero 500 ptos
Ts = 1/Fs; 
x=0:Ts:(N-1)*Ts; %eje discreto
f0=10; %hz
f1=0.1;%hz
y=[sin(2*pi*f0*x(1:floor(N/2))),sin(2*pi*f1*x(floor(N/2):end))];
% Mexican Hat escalable
%mex_hat = @(t, a, b) (1 - ((t - b)/a).^2) .* exp(-((t - b).^2) / (2*a^2)) / sqrt(abs(a));

% Morlet real escalable
%morlet = @(t, a, b) cos(5*((t - b)/a)) .* exp(-((t - b).^2)/(2*a^2)) / sqrt(abs(a));

% STFT (por ventana de gabor)
%spectrogram(y)

% CWT (wavelet continua)
% figure(1)
% scales = 0.5:1:100; % usar solo escalas del 1 al 50
% coef = cwt(y, scales, 'morl'); 
% imagesc(coef); axis xy; colorbar



plot_dwt_levels(y, Fs, 12, 'db4');
