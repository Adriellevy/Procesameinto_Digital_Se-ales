%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%% Ejercicio Nro. 1) de la Guía_5
%%%
%%% Unidad N° 5
%%%
%%% Análisis de Señales con ventanas de Gabor
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

Ejercicio_1(f,fs,al,Nd);   % Llamado a la Transformada de Gabor


