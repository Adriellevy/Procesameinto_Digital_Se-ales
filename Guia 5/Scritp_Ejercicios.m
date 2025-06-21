% % 1) Test Algoritmo Transformada de Gabor
% N=512;
% T0=2;
% fs=N/T0;  % 256 Hz en este caso
% Ts=1/fs;
% al=1/512;  % DispersiÃ³n de la ventana de Gabor
% t1=0:Ts:T0/2-Ts;
% t2=T0/2:Ts:T0-Ts;
% f1=2*sin(2*pi*64*t1);
% f2=sin(2*pi*16*t2);
% t=[t1 t2];
% f=[f1 f2];  % tonos de 16 y 64 Hz a diferentes tiempos 
% Nd=64; % Cantidad de desplazamientos de la ventana de Gabor
% otra forma de hacerlo 
Ts = 1/256;             % Periodo de muestreo (256 Hz como el segundo ejemplo)
x = 0:Ts:2-Ts;          % Vector de tiempo total de 2 segundos
mid = floor(length(x)/2);

y1 = sin(2 * pi * 64 * x(1:mid));
y2 = sin(2 * pi * 16 * x(mid+1:end));
y = [y1 y2];

[wt,tnew,fnew] = TransformadaGabor(f,fs,1/512,Nd);

figure;
subplot(2,1,1);
plot(t,f);
xlabel('t (s)');
ylabel('y(t)');
title('Senial original');
grid on;

subplot(2,1,2);
mesh(fnew, tnew, wt);
axis([min(fnew) max(fnew) min(tnew) max(tnew) min(min(wt)) max(max(wt))]);
xlabel('f (Hz)');
ylabel('t (s)');
zlabel('Amplitud');
title('Transformada de Gabor');
rotate3d;

%---------------------------------------------------------------
% % 2) Señal con 3 senos para probar la transformada de gabor
% % Parámetros
% N = 512;
% fs = 128;
% t = linspace(0, 4, N);
% x = zeros(size(t));
% 
% % Construcción de la señal por tramos
% for i = 1:length(t)
%     if t(i) < 1
%         x(i) = sin(2*pi*32*t(i));
%     elseif t(i) < 2
%         x(i) = sin(2*pi*16*t(i));
%     else
%         x(i) = sin(3*pi*8*t(i));
%     end
% end
% 
% % Aplicar Transformada de Gabor
% Alfa = 1/256;
% ModulacionFrecuencial = 100; % más modulaciones para ver evolución temporal
% [G, t_new, f_new] = TransformadaGabor(x, fs, Alfa, ModulacionFrecuencial);
% 
% % Graficar
% figure;
% imagesc(t_new, f_new, G');
% axis xy;
% xlabel('Tiempo (s)');
% ylabel('Frecuencia (Hz)');
% title('Transformada de Gabor - Señal por Tramos');
% colorbar;
%B) Por ventana rectangular 
% Parámetros
% N = 1024;
% fs = 512;
% t = linspace(0, N/fs, N);
% x = double(t >= 1 & t <= 1.1);  % señal rectangular
% 
% % Aplicar Transformada de Gabor
% Alfa = 641;
% ModulacionFrecuencial = 100;
% [G, t_new, f_new] = TransformadaGabor(x, fs, Alfa, ModulacionFrecuencial);
% 
% % Graficar
% figure;
% imagesc(t_new, f_new, G');
% axis xy;
% xlabel('Tiempo (s)');
% ylabel('Frecuencia (Hz)');
% title('Transformada de Gabor - Señal Rectangular');
% colorbar;

%---------------------------------------------------------------
% 3) Demostracion que es lo mismo señal con menos pasos
% a = 2;  % escala
% 
% % Opción A: g(at), intervalo [-5a, 5a], muestreo en enteros
% tA = -5*a:1:5*a;
% gA = exp(-tA.^2);  % ventana gaussiana
% 
% % Opción B: g(t), intervalo [-5,5], muestreo con paso 1/a
% tB = -5:1/a:5;
% gB = exp(-tB.^2);  % misma ventana sin escalar
% 
% % Graficar
% figure;
% subplot(2,1,1);
% stem(tA, gA);
% title('a) g(at), t ? [-5a, 5a], muestreo entero');
% 
% subplot(2,1,2);
% stem(tB, gB);
% title('b) g(t), t ? [-5, 5], muestreo 1/a');
%---------------------------------------------------------------
% 4) Onditas Mex y morlet
% N_values = [81, 65, 49, 33, 17];
% 
% figure;
% for i = 1:length(N_values)
%     N = N_values(i);
%     t = linspace(-8, 8, N);
% 
%     % Ondita Sombrero Mexicano
%     psi_mex = (1 - t.^2) .* exp(-t.^2 / 2);
%     
%     % Ondita Morlet
%     psi_morlet = cos(5*t) .* exp(-t.^2 / 2);
% 
%     % FFTs
%     ang_mex = abs(fftshift(fft(psi_mex)));
%     ang_mor = abs(fftshift(fft(psi_morlet)));
% 
%     % Plot
%     subplot(2, length(N_values), i);
%     plot(ang_mex);
%     title(['Mex N=' num2str(N)]);
%     
%     subplot(2, length(N_values), i + length(N_values));
%     plot(ang_mor);
%     title(['Morlet N=' num2str(N)]);
% end


%---------------------------------------------------------------
% 5) Señal escalonada que define la ondícula Haar
% t = linspace(0, 1, 1024);
% psi = double(t >= 0) - 2*double(t >= 0.5) + double(t >= 1);
% 
% figure;
% plot(t, psi);
% title('Ondita Madre de Haar (Interpolada N = 1024)');
% xlabel('t');
% ylabel('\psi(t)');

%---------------------------------------------------------------
% 6) Filtro Daubechies 4
% 
% % Define el nombre de la wavelet
% wname = 'db4';
% 
% % Número de puntos para la representación
% N = 1024;
% 
% % Obtener la función wavelet (ondita) y escalamiento
% [phi, psi, xval] = wavefun(wname, 10);  % 10 niveles de descomposición
% 
% % Graficar la ondita (función wavelet madre)
% figure;
% plot(xval, psi);
% title('Onda madre de Daubechies 4 (db4)');
% xlabel('Tiempo');
% ylabel('\psi(t)');
% grid on;

% % Wavedec Descomosicion de las frecuencias en mitades. 
% 
% % Limpieza
% clear; close all; clc;
% 
% % Parámetros de la señal
% Fs = 1000;           % Frecuencia de muestreo (Hz)
% t = 0:1/Fs:1-1/Fs;   % Vector de tiempo (1 segundo)
% 
% % Señal compuesta: suma de senos de diferentes frecuencias
% s1 = sin(2*pi*10*t);    % 10 Hz
% s2 = sin(2*pi*50*t);    % 50 Hz
% s3 = sin(2*pi*200*t);   % 200 Hz
% 
% % Señal total
% x = s1 + s2 + s3;
% 
% % Mostrar la señal original
% figure;
% plot(t, x);
% title('Señal original: suma de senos');
% xlabel('Tiempo (s)');
% ylabel('Amplitud');

% 


%---------------------------------------------------------------
% % Descomposición wavelet
% nivel = 5;                         % Número de niveles de descomposición
% wavelet = 'db4';                   % Tipo de wavelet (Daubechies 4)
% [C, L] = wavedec(x, nivel, wavelet);
% 
% % Extraer los coeficientes de detalle y aproximación
% A5 = appcoef(C, L, wavelet, 5);    % Aproximación en nivel 5
% D5 = detcoef(C, L, 5);             % Detalle nivel 5
% D4 = detcoef(C, L, 4);             % Detalle nivel 4
% D3 = detcoef(C, L, 3);             % Detalle nivel 3
% D2 = detcoef(C, L, 2);             % Detalle nivel 2
% D1 = detcoef(C, L, 1);             % Detalle nivel 1
% 
% % Mostrar los coeficientes
% figure;
% subplot(6,1,1); plot(A5); title('Aproximación nivel 5');
% subplot(6,1,2); plot(D5); title('Detalle nivel 5');
% subplot(6,1,3); plot(D4); title('Detalle nivel 4');
% subplot(6,1,4); plot(D3); title('Detalle nivel 3');
% subplot(6,1,5); plot(D2); title('Detalle nivel 2');
% subplot(6,1,6); plot(D1); title('Detalle nivel 1');




