%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Ejercicio 9_a de Cálculo de CWT con pulso rectangular y ondita Haar
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



function X=Ejercicio_9_a_CWT

x=[8 8 8 8 0 0 0 0];
n=0:length(x)-1;

H=figure(1);
set(H,'NumberTitle','off','Menubar','none','name','Ejercicio 10. Ondita Continua. Pulso de energía','position',[60 30 1000 700]);
subplot(121)
stem(n,x)
xlabel('n');
ylabel('x[n]');
title('Pulso de energía')
subplot(122)
Esc=1:16;
X=cwt(x,Esc,'haar','plot');
Y=abs(round(X)); % Vista de la matriz como en el libro para graficar
M=[Esc' Y]