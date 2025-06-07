%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Ejeercicio 9b de Cálculo de CWT con 2 tonos senoidales puros de energía
%%% utilizando Ondita Morlet
%%%
%%% Unidad N° 5
%%%
%%% Procesamiento Digital de Señales
%%%
%%% Dr. Ing. Franco Martin Pessana
%%% 
%%% Facultad de Ingeniería y Ciencias Exactas y Naturales
%%% Universidad Favaloro
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function X=Ejercicio_9_b_CWT

fs=500;
Ts=1/fs;
N=100;
t=0:Ts:(N-1)*Ts;
f01=200;
f02=50;
x1=3*sin(2*pi*f01*t(1:floor(N/2)));
x2=sin(2*pi*f02*t(floor(N/2):N-1));
x=[x1 x2];
H1=figure(1);
set(H1,'NumberTitle','off','Menubar','none','name','Ejemplo #2. Ondita Continua. Pulsos de energía Senoidales','position',[60 30 1000 700]);
subplot(121)
plot(t,x);
xlabel('t');
ylabel('x(t)');
title('pulsos de senoides')
grid;
subplot(122)
Esc=1:32    ;
X=cwt(x,Esc,'morl','plot');
Y=abs(round(X));
H2=figure(2);
set(H2,'NumberTitle','off','Menubar','none','name','Ejemplo #2. Vista 3D','position',[100 100 600 600]);
surf(abs(X));
xlabel('n');
ylabel('esc');
title('|X(b,a)|')
rotate3d;
