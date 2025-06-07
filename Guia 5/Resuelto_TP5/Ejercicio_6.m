%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Ejercicio 6) TP_5
%%% 
%%% Construcción de Ondita DB4 a partir de h'[n] y l'[n]
%%%
%%% Procesamiento Digital de Señales
%%%
%%% Dr. Ing. Franco Martin Pessana
%%% 
%%% Facultad de Ingeniería y Ciencias Exactas y Naturales
%%% Universidad Favaloro
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ejercicio_6

lp=[.4830 .8365 .2241 -.1294];
n=0:length(lp)-1;
hp=(-1).^n.*lp(length(lp)-n);
%hp=[-.1294 -.2241 .8365 -.4830];

Level=8;
H1=figure(1);
set(H1,'NumberTitle','off','Menubar','none','name','Ondita Discreta Daubechies 4','position',[60 30 1000 700]);
Ondita_DB4=hp;
[W_Ondita_DB4,w]=freqz(Ondita_DB4);
subplot(121)
stem(n,Ondita_DB4)
ylabel('phi[n]')
title('Ondita DB4');
grid
subplot(122)
plot(w,abs(W_Ondita_DB4));
xlabel('w (rad)');
ylabel('|W_DB4|')
title('Respuesta Frecuencial')
grid
disp('Presione cualquier tecla para continuar con el escalamiento temporal...')
pause
for k=2:Level
    Ondita_DB4=upsample(Ondita_DB4,2);
    Ondita_DB4=Ondita_DB4(1:length(Ondita_DB4)-1); % Ya que upsample agrega un cero al final
    Ondita_DB4=conv(Ondita_DB4,lp);
    [W_Ondita_DB4,w]=freqz(Ondita_DB4);
    %length(Ondita_DB4)
    subplot(121)
    n=0:length(Ondita_DB4)-1;
    stem(n,Ondita_DB4)
    ylabel('phi[n]')
    title('Ondita DB4');
    grid
    subplot(122)
    plot(w,abs(W_Ondita_DB4));
    xlabel('w (rad)');
    ylabel('|W_DB4|')
    title('Respuesta Frecuencial')
    grid
    disp('Presione cualquier tecla para continuar con el escalamiento temporal...')
    pause
end
subplot(121)
hold on
plot(n,Ondita_DB4,'r*')