%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%% Ejercicio Nro. 10) de la Guía_5
%%% Detalles y aproximaciones de una señal digital de 16 bits
%%% (cada bit compuesto por 8 valores)
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

si=[1 0 1 1 0 0 1 0 1 1 1 0 1 0 0 1]; % Señal de 16 bits

% Construcción de la señal digital

Byte=8;
N=length(si);
Level=4;             % Cantidad de niveles de DWT
Si=zeros(1,Byte*N);
n=0:length(Si)-1;   % Variable independiente de entrada
for k=1:N
    if si(k)
        Si(Byte*(k-1)+1:Byte*k)=ones(1,8);
    end
end
H1=figure(1);
set(H1,'position',[300 200 800 500],'Menubar','none',...
        'NumberTitle','off','name',...
        'Pulso binario de 16 bits y 8 bits por byte');
stem(n,Si)
% save senial_binaria.mat Si;
H2=figure(2);
set(H2,'position',[20 30 1200 700],'Menubar','none',...
        'NumberTitle','off','name',...
        'Ejemplo de Denoising con DWT. PDS. FICEN');
subplot(5,2,1)
stem(n,Si)
subplot(5,2,2)
Si_CWT=cwt(Si,1:1:10,'haar','plot');
[Si_C_DWT,L]=wavedec(Si,Level,'haar');
App_Level=wrcoef('a',Si_C_DWT,L,'haar',Level);
Det_k=zeros(Level,length(Si));
for k=1:Level
    Det_k(k,:)=wrcoef('d',Si_C_DWT,L,'haar',k);
end
subplot(5,2,4)
stem(n,Det_k(1,:))
ylabel('Det_1')
subplot(5,2,6)
stem(n,Det_k(2,:))
ylabel('Det_2')
subplot(5,2,8)
stem(n,Det_k(3,:))
ylabel('Det_3')
subplot(5,2,9)
stem(n,App_Level)
ylabel('App_4')
xlabel('Aproximaciones')
subplot(5,2,10)
stem(n,Det_k(Level,:))
ylabel('Det_4')
xlabel('Detalles')

subplot(5,2,7)
App_3=App_Level+Det_k(4,:);
stem(n,App_3)
ylabel('App_3')
subplot(5,2,5)
App_2=App_3+Det_k(3,:);
stem(n,App_2)
ylabel('App2')
subplot(5,2,3)
App_1=App_2+Det_k(2,:);
stem(n,App_1)
ylabel('App_1')