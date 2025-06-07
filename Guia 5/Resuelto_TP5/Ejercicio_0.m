%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%% Ejercicio Nro. 0) de la Guía_5
%%%
%%% Unidad N° 5
%%%
%%% Análisis de Señales con Transformada de Fourier
%%% 
%%% Señal de Prueba (2 tonos senoidales) para Gabor
%%%
%%% Procesamiento Digital de Señales
%%%
%%% Dr. Ing. Franco Martin Pessana
%%% 
%%% Facultad de Ingeniería y Ciencias Exactas y Naturales
%%% Universidad Favaloro
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
T0=2;
Fs=256;  % 256 Hz en este caso
Ts=1/Fs;
N=T0/Ts;
t1=0:Ts:T0/2-Ts;
t2=T0/2:Ts:T0-Ts;
f1=2*sin(2*pi*64*t1);
f2=sin(2*pi*16*t2);
t=[t1 t2];
f=[f1 f2];  % tonos de 16 y 64 Hz a diferentes tiempos 
fd=0:1/T0:Fs-1/T0;   % Frecuencia en Hz
F=fft(f);

H1=figure(1);
set(H1,'position',[20 100 1200 650],'Menubar','none',...
        'NumberTitle','off','name','TP#5. Transformada Discreta de Fourier de pulsos senoidales');
subplot(211)
plot(t,f)
xlabel('t(s)');
ylabel('f(t)');
grid
subplot(212)
plot(fd,abs(F))
axis([min(fd) max(fd) min(abs(F)) max(abs(F))])
xlabel('f (Hz)');
ylabel('abs(F)');
grid


