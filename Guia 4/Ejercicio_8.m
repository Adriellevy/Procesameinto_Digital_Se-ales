%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%% Ejercicio Nro. 8) de la Guía_4
%%% 
%%% Unidad N° 4
%%%                                                            
%%% Cancelación Adaptativa de ruido en Resonancia Magnética. 
%%% Se utiliza el modelo predictivo que toma como señal deseada el ruido
%%% mas la voz del paciente y como entrada el ruido del resonador.
%%% Luego se la resta la salida del filtro adaptatico a la señal deseada.    
%%%  
%%% Sintaxis:
%%%                
%%%            [B]=Ejercicio_8(L)
%%%
%%%
%%% Procesamiento Digital de Señales
%%%
%%% Dr. Ing. Franco Martin Pessana
%%% 
%%% Facultad de Ingeniería y Ciencias Exactas y Naturales
%%% Universidad Favaloro
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B,A]=Ejercicio_8(L)

if nargin<1
   L=500; % El orden óptimo es XXX
end

load 'fMRI_mic.mat';
xk=mic2; % Ruido del resonador
dk=mic1; % Ruido del resonador mas señal de voz
N=length(xk);
Fs=fs;
Ts=1/Fs;			% Período de muestreo
t=0:Ts:(N-1)*Ts;	% Vector de tiempos

[B,A,yk,ek]=ARMA_Adaptativo(xk,dk,L,0); %% Llamada a IIR Adaptativo como FIR Q=0

H=figure(1);
set(H,'NumberTitle','off','Menubar','none','name','FIR Adaptativo Ejercicio 8','position',[20 30 1000 700]);
subplot(311)
plot(t,dk,'b',t,xk,'r');
ylabel('xk red; dk blue');
title('Ruido en Resonador y voz sumergida en ruido MRI')
subplot(312)
plot(t,yk);
ylabel('yk (salida Filtro)');
title('Salida del Filtro FIR Adaptativo en modalidad Cancelador')
subplot(313)
plot(t,ek);
ylabel('ek (Voz sin ruido)');
xlabel('t [s]')
title('Cancelación de ruido de resonador en voz del paciente')

disp('%%%%%%%%%% Voz del Paciente con Ruido de Resonador (80 segundos aprox): %%%%%%%%%%%')
sound(dk(1:10000),Fs);
disp('%%%%%%%% Presione ENTER para continuar...%%%%%%%%%%%%')
pause
disp('%%%%%%%%%% Voz del Paciente sin Ruido de Resonador (80x` segundos aprox): %%%%%%%%%%%')
sound(ek(1:30000),Fs);

H2=figure(2);
set(H2,'NumberTitle','off','Menubar','none','name','FIR Adaptativo Ejercicio 8',...
    'position',[60 70 1000 700]);
Ek=fft(ek);
Af=(Fs/N);
f=0:Af:Fs-Af;
f1=0.05e4;
f2=0.09e4;
N1=round(N*f1/Fs);
N2=round(N*f2/Fs);
H=ones(size(Ek));
H(N1:N2)=0;
H(N-N2:N-N1)=0;
Ek_filt=Ek.*H;
ek_filt=real(ifft(Ek_filt));
subplot(2,1,1)
plot(f,abs(fft(xk)));
subplot(2,1,2)
%plot(f,abs(Ek_filt));
plot(f,abs(Ek),'r',f,max(abs(Ek))*H,'b');
%sound(ek_filt,Fs)
end