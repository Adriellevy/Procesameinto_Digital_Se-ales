%%%%%%%%%%%%%%%%%%5555%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%% Ejercicio Nro. 5)
%%%
%%% Unidad N° 5. Transformada Ondita
%%%
%%% Construcción de Ondita Haar a partir de h'[n] y l'[n]
%%%
%%% Procesamiento Digital de Señales
%%%
%%% Dr. Ing. Franco Martin Pessana
%%% 
%%% Facultad de Ingeniería y Ciencias Exactas y Naturales
%%% Universidad Favaloro
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ejercicio_5

hp=[1 -1];
lp=[1 1];

Level=10;
H1=figure(1);
set(H1,'NumberTitle','off','Menubar','none','name','Ondita Continua Haar','position',[60 30 1000 700]);
Ondita_Haar=hp;
[W_Ondita_Haar,w]=freqz(Ondita_Haar);
subplot(121)
n=0:length(Ondita_Haar)-1;
stem(n,Ondita_Haar)
ylabel('phi[n]')
title('Ondita Haar');
grid
subplot(122)
plot(w,abs(W_Ondita_Haar));
xlabel('w (rad)');
ylabel('|W_Haar|')
title('Respuesta Frecuencial')
grid
disp('Presione cualquier tecla para continuar con el escalamiento temporal...')
pause
for k=2:Level
    Long=2^k;
    Ondita_Haar=wkeep(conv(dyadup(Ondita_Haar),lp),Long);
    [W_Ondita_Haar,w]=freqz(Ondita_Haar);
    subplot(121)
    n=0:length(Ondita_Haar)-1;
    stem(n,Ondita_Haar)
    ylabel('phi[n]')
    title('Ondita Haar');
    grid
    subplot(122)
    plot(w,abs(W_Ondita_Haar)/sqrt(2^k));
    xlabel('w (rad)');
    ylabel('|W_Haar|')
    title('Respuesta Frecuencial')
    grid
    disp('Presione cualquier tecla para continuar con el escalamiento temporal...')
    pause
end
subplot(121)
hold on
plot(n,Ondita_Haar,'r*')