%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%% Ejercicio Nro. 2a)
%%%
%%% Unidad N° 5. Transformada Ondita
%%%
%%% Análisis de Señales con ventanas de Gabor
%%% 
%%% Señal de Prueba (3 tonos senoidales) para Gabor
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
N=512;
T0=4;
fs=N/T0;  % 256 Hz en este caso
Ts=1/fs;
al=1/256;  % Dispersión de la ventana de Gabor
t1=0:Ts:1-Ts;
t2=1:Ts:2-Ts;
t3=2:Ts:T0-Ts;

f1=3*sin(2*pi*8*t1);
f2=2*sin(2*pi*16*t2);
f3=sin(2*pi*32*t3);
t=[t1 t2 t3];
f=[f1 f2 f3];  % tonos de 8, 16 y 32 Hz a diferentes tiempos 
Nd=64; % Cantidad de desplazamientos de la ventana de Gabor

Ejercicio_1(f,fs,al,Nd);   % Llamado a la Transformada de Gabor