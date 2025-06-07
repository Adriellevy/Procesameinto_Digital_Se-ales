%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Ejemplo de Cálculo de DWT chirp signal (senoide con incremento de
%%% frecuencia, utilizando hasta 4 nivel de la DWT utilizando 
%%% Ondita Daubechies 10
%%%
%%% MSc Ing. Franco Martin Pessana
%%% Procesamiento Avanzado de Señales
%%% Facultad Regional Buenos Aires
%%% Universidad Tecnológica Nacional
%%% Doctorado en Ingeniería
%%% Ciclo Lectivo 2012
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Si_C_DWT=Ejemplo_4_DWT()

Level=4;

fs=1000;
Ts=1/fs;
Tf=2;
Fi=0;
Ff=150;
t=0:Ts:Tf-Ts;                    
Si=chirp(t,Fi,Tf,Ff);

H1=figure(1);
set(H1,'position',[300 200 800 500],'Menubar','none',...
        'NumberTitle','off','name',...
        'Chirp Signal de 2 segundos de duración');
plot(t,Si)
H2=figure(2);
set(H2,'position',[20 30 1200 700],'Menubar','none',...
        'NumberTitle','off','name',...
        'Ejemplo#4. Aproximaciones y detalles de 4 nivel de DWT de Chirp Signal.');
subplot(5,2,1)
plot(t,Si)
subplot(5,2,2)
Si_CWT=cwt(Si,1:1:64,'morl','plot');
[Si_C_DWT,L]=wavedec(Si,Level,'db10');
App_Level=wrcoef('a',Si_C_DWT,L,'db10',Level);
Det_k=zeros(Level,length(Si));
for k=1:Level
    Det_k(k,:)=wrcoef('d',Si_C_DWT,L,'db10',k);
end
subplot(5,2,4)
plot(t,Det_k(1,:))
ylabel('Det_1')
subplot(5,2,6)
plot(t,Det_k(2,:))
ylabel('Det_2')
subplot(5,2,8)
plot(t,Det_k(3,:))
ylabel('Det_3')
subplot(5,2,9)
plot(t,App_Level)
ylabel('App_4')
xlabel('Aproximaciones')
subplot(5,2,10)
plot(t,Det_k(Level,:))
ylabel('Det_4')
xlabel('Detalles')
subplot(5,2,7)
App_3=App_Level+Det_k(4,:);
plot(t,App_3)
ylabel('App_3')
subplot(5,2,5)
App_2=App_3+Det_k(3,:);
plot(t,App_2)
ylabel('App2')  
subplot(5,2,3)
App_1=App_2+Det_k(2,:);
plot(t,App_1)
ylabel('App_1')