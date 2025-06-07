%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%% Ejercicio Nro. 4)
%%%
%%% Unidad N° 5. Transformada Ondita
%%%
%%% Construcción de Ondita Sombrero Mexicano (Ondita Cruda)
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

function Ejercicio_4_Ondita_Mexh

t=-8:1:8;
Ondita_Mexh=2/(sqrt(3)*pi^0.25).* exp(-t.^2/2).*(1-t.^2);
Level=8;
H1=figure(1);
set(H1,'NumberTitle','off','Menubar','none','name','Ondita Continua Sombrero Mexicano','position',[60 30 1000 700]);
[W_Ondita_Mexh,w]=freqz(Ondita_Mexh);
subplot(121)
stem(t,Ondita_Mexh)
axis([min(t) max(t) min(Ondita_Mexh) max(Ondita_Mexh)]);
ylabel('phi[n]')
title('Ondita Sombrero Mexicano');
grid
subplot(122)
plot(w,abs(W_Ondita_Mexh));
xlabel('w (rad)');
ylabel('|W_Mexh|')
title('Respuesta Frecuencial')
grid
disp('Presione cualquier tecla para continuar con el escalamiento temporal...')
pause
for k=2:Level;
    t=-8:1/k:8;
    Ondita_Mexh=2/(sqrt(3)*pi^0.25).* exp(-t.^2/2).*(1-t.^2);
    [W_Ondita_Mexh,w]=freqz(Ondita_Mexh);
    subplot(121)
    stem(t,Ondita_Mexh)
    axis([min(t) max(t) min(Ondita_Mexh) max(Ondita_Mexh)]);
    ylabel('phi[n]')
    title('Ondita Sombrero Mexicano');
    grid
    subplot(122)
    plot(w,abs(W_Ondita_Mexh)/sqrt(k));
    xlabel('w (rad)');
    ylabel('|W_Mexh|')
    title('Respuesta Frecuencial')
    grid
    disp('Presione cualquier tecla para continuar con el escalamiento temporal...')
    pause
end