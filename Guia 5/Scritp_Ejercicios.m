
N=512;
T0=2;
fs=N/T0;  % 256 Hz en este caso
Ts=1/fs;
al=1/512;  % Dispersión de la ventana de Gabor
t1=0:Ts:T0/2-Ts;
t2=T0/2:Ts:T0-Ts;
f1=2*sin(2*pi*64*t1);
f2=sin(2*pi*16*t2);
t=[t1 t2];
f=[f1 f2];  % tonos de 16 y 64 Hz a diferentes tiempos 
Nd=64; % Cantidad de desplazamientos de la ventana de Gabor
% otra forma de hacerlo 
% Ts = 1/256;             % Periodo de muestreo (256 Hz como el segundo ejemplo)
% x = 0:Ts:2-Ts;          % Vector de tiempo total de 2 segundos
% mid = floor(length(x)/2);
% 
% y1 = sin(2 * pi * 64 * x(1:mid));
% y2 = sin(2 * pi * 16 * x(mid+1:end));
% y = [y1 y2];

[wt,tnew,fnew] = TransformadaGabor(f,fs,1/512,Nd);

figure;
subplot(2,1,1);
plot(t,f);
xlabel('t (s)');
ylabel('y(t)');
title('Señal original');
grid on;

subplot(2,1,2);
mesh(fnew, tnew, wt);
axis([min(fnew) max(fnew) min(tnew) max(tnew) min(min(wt)) max(max(wt))]);
xlabel('f (Hz)');
ylabel('t (s)');
zlabel('Amplitud');
title('Transformada de Gabor');
rotate3d;
