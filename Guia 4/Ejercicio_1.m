%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%% Ejercicio Nro. 1) de la Guía_4
%%%
%%% Unidad N° 4
%%%
%%% Filtro FIR adaptativo
%%%
%%% Procesamiento Digital de Señales
%%%
%%% Dr. Ing. Franco Martin Pessana
%%% 
%%% Facultad de Ingeniería y Ciencias Exactas y Naturales
%%% Universidad Favaloro
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

Nro=1000;
N=5;
k=(0:Nro)';
xk=sin(2*pi*k/N);			% Señal de Entrada
dk=2*cos(2*pi*k/N);			% Señal deseada
P=1;                        % Orden del Filtro FIR
% [Bk,yk,ek]=FIR_Adaptativo(xk,dk,P);
[Bk,Ak,yk,ek]=ARMA_Adaptativo(xk,dk,P,0);
Bk
Ak
BKT=[2*cos(2*pi/N)/sin(2*pi/N) -2/sin(2*pi/N)]'
H1=figure(1);
set(H1,'position',[20 30 1000 700],'Menubar','none',...
        'NumberTitle','off','name','Ejercicio 1. Filtro FIR Adaptativo');
subplot(311)
plot(k,xk);
ylabel('xk');
subplot(312)
plot(k,yk);
ylabel('yk');
subplot(313)
plot(k,ek);
ylabel('ek');
xlabel('k')
H2=figure(2);
set(H2,'position',[60 30 1000 700],'Menubar','none',...
        'NumberTitle','off','name','Superficie de Minimización');
w0=-10:0.2:10;
Nx=size(w0,2);
B0=ones(Nx,1)*w0;% Matriz con variaciones de wo
B1=B0';				% Matriz con variaciones de w1
Ek=2+0.5*(B0.^2+B1.^2)+B0.*B1*cos(2*pi/N)+2*B1*sin(2*pi/N);
mesh(B0,B1,Ek)
xlabel('b0')
ylabel('b1')
zlabel('E[ek*ek]')
rotate3d;
H3=figure(3);
set(H3,'position',[100 30 1000 700],'Menubar','none',...
        'NumberTitle','off','name','Curvas de Nivel del Error');
contour(B0,B1,Ek,30)
xlabel('b0')
ylabel('b1')