 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%% Ejercicio Nro. 4)
%%%
%%% Unidad N° 5. Transformada Ondita
%%%
%%% Construcción de Ondita Morlet (Ondita Cruda)
%%%
%%% MSc Ing. Franco Martin Pessana
%%% Procesamiento Avanzado de Señales
%%%
%%% Dr. Ing. Franco Martin Pessana
%%% 
%%% Facultad de Ingeniería y Ciencias Exactas y Naturales
%%% Universidad Favaloro
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ejercicio_4_Ondita_Morlet

t=-8:1:8;
Ondita_Morl=exp(-t.^2/2).*cos(5*t);

Level=8;
H1=figure(1);
set(H1,'NumberTitle','off','Menubar','none','name','Ondita Continua Morlet','position',[60 30 1000 700]);
[W_Ondita_Morl,w]=freqz(Ondita_Morl);
subplot(121)
stem(t,Ondita_Morl)
axis([min(t) max(t) min(Ondita_Morl) max(Ondita_Morl)]);
ylabel('phi[n]')
title('Ondita Morlet');
grid
subplot(122)
plot(w,abs(W_Ondita_Morl));
xlabel('w (rad)');
ylabel('|W_Morl|')
title('Respuesta Frecuencial')
grid
disp('Presione cualquier tecla para continuar con el escalamiento temporal...')
pause
for k=2:Level;
    t=-8:1/k:8;
    Ondita_Morl=exp(-t.^2/2).*cos(5*t);
    [W_Ondita_Morl,w]=freqz(Ondita_Morl);
    subplot(121)
    stem(t,Ondita_Morl)
    axis([min(t) max(t) min(Ondita_Morl) max(Ondita_Morl)]);
    ylabel('phi[n]')
    title('Ondita Morlet');
    grid
    subplot(122)
    plot(w,abs(W_Ondita_Morl)/sqrt(k));
    xlabel('w (rad)');
    ylabel('|W_Morl|')
    title('Respuesta Frecuencial')
    grid
    disp('Presione cualquier tecla para continuar con el escalamiento temporal...')
    pause
end