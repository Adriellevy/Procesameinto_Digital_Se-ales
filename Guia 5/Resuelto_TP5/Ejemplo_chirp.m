%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Ejercicio 7) TP_5
%%% 
%%% Evaluación de Ondita Sombrero Mexicano y Morlet con pulsos senoidales
%%%
%%% Procesamiento Digital de Señales
%%%
%%% Dr. Ing. Franco Martin Pessana
%%% 
%%% Facultad de Ingeniería y Ciencias Exactas y Naturales
%%% Universidad Favaloro
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c=Ejemplo_chirp

fs=1e3;
Ts=1/fs;
N=1024;
x=zeros(1,N);
t=0:Ts:(N-1)*Ts;

x=chirp(t);

H1=figure(1);
set(H1,'NumberTitle','off','Menubar','none','name','Pulsos Senoidales con Ondita Sombrero Mexicano','position',[60 30 1000 700]);
subplot(121)
plot(t,x);
axis([min(t) max(t) min(x) max(x)])
xlabel('t (s)');
ylabel('x(t)')
title('Señal temporal')
grid
subplot(122)
c=cwt(x,1:12,'mexh','plot');
%axis([min(t) max(t) min(x) max(x)])
xlabel('t (s)');
ylabel('escala')
title('Magnitud de la CWT')
grid

H2=figure(2);
set(H2,'NumberTitle','off','Menubar','none','name','Pulsos Senoidales con Ondita Morlet','position',[100 10 1000 700]);
subplot(121)
plot(t,x);
axis([min(t) max(t) min(x) max(x)])
xlabel('t (s)');
ylabel('x(t)')
title('Señal temporal')
grid
subplot(122)
c=cwt(x,1:22,'morl','plot');
%axis([min(t) max(t) min(x) max(x)])
xlabel('t (s)');
ylabel('escala')
title('Magnitud de la CWT')
grid