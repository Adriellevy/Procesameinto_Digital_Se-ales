%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%% Ejercicio 11 TP#2. Filtro Notch Digital
%%%
%%% Filtro Notch
%%%
%%% Procesamiento Digital de Señales
%%%
%%% Dr. Ing. Franco Martin Pessana
%%% 
%%% Facultad de Ingeniería y Ciencias Exactas y Naturales
%%% Universidad Favaloro
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cómputo de Filtro Notch

clear
clc

k = 512;  %puntos para graficado de espectros
fo=50;  % Frecuencia de línea a eliminar
Fs=1000;  % Trecuencia de muestreo
wo=2*pi*fo/Fs; % Ángulo de ubicación de los ceros del filtro IIR
r = 0.9; % Magnitud de los ceros
num = conv([1 -exp(-1j*wo)],[1 -exp(1j*wo)])
den = conv([1 -r*exp(-1i*wo)],[1 -r*exp(1i*wo)])
[h,w] = freqz(num, den, k);
H=figure(1);
set(H,'NumberTitle','off','Menubar','none','name',...
    'Filtro Notch. Ejercicio#11. FICEN UF','position',[20 50 1400 750]);
% Graficado de la Respuesta en Frecuencia
subplot(2,2,1)
plot(w,real(h));grid
title('parte real')
xlabel('omega [rad]'); ylabel('Amplitud')
subplot(2,2,2)
plot(w,imag(h));grid
title('Parte imaginaria')
xlabel('omega [rad]'); ylabel('Amplitud')
subplot(2,2,3)
plot(w,abs(h));grid
title('Espectro de Magnitud')
xlabel('omega [rad]'); ylabel('Magnitud')
subplot(2,2,4)
plot(w,angle(h));grid
title('espectro de Fase')
xlabel('omega [rad]'); ylabel('Fase, radianes')
N=150;
Ts=1/Fs;
t=0:Ts:(N-1)*Ts;
x=5*sin(2*pi*fo*t);
y=filter(num,den,x);
H=figure(2);
set(H,'NumberTitle','off','Menubar','none','name',...
    'Entrada y Salida del Filtro Notch. Ejercicio#11. FICEN UF',...
    'position',[100 20 1400 750]);
subplot(211)
plot(t,x)
grid
xlabel('t [s]')
ylabel('x(t)')
subplot(212)
plot(t,y)
xlabel('t [s]')
ylabel('y(t)')
grid
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cálculo te tiempo de establecimiento del Filtro Notch
Max=max(y);
Umb=(abs(y)<0.01*Max)*Max;
for n=N:-1:1
    if Umb(n)==0
        IndTe = n;
        break
    end   
end
Test=t(IndTe);    % Tiempo de Establecimiento
plot(t,Umb,'g');
disp('Tiempo de establecimiento (s) =')
Test