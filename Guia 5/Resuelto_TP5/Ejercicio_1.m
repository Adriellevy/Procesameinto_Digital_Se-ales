%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%% Ejercicio Nro. 1)
%%%
%%% Unidad N° 5. Transformada Ondita
%%%
%%% Análisis de Señales con ventanas de Gabor. Visualizado Tridimensional
%%%
%%%           G_F_GW=Ejercicio_1(f,Fs,al,Nd)
%%%
%%% Parámetros de entrada:
%%%             Señal temporal: f
%%%             Frecuencia de muestreo: Fs
%%%             Dispersion de ventana de Gabor: al
%%%             Cantidad de Desplazamientos temporales: Nd
%%%
%%% Parámetros de salida:
%%%             magnitud de la transformada de Gabor: G_F_GW
%%%
%%% Procesamiento Digital de Señales
%%%
%%% Dr. Ing. Franco Martin Pessana
%%% 
%%% Facultad de Ingeniería y Ciencias Exactas y Naturales
%%% Universidad Favaloro
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G_F_GW=Ejercicio_1(f,Fs,al,Nd)

N=length(f);            % Longitud de la señal f
Ts=1/Fs;                % Período de muestreo de la señal f
t=0:Ts:(N-1)*Ts;        % eje temporal discreto
T0=N*Ts;                % Ventana de observación temporal
At=T0/Nd;   % Intervalos de tiempo para corrimiento de la Ventana de Gabor

G_F_GW=zeros(Nd+1,length(t));  % Inicialización de Matriz de Gabor (Modulo)
for k=1:Nd+1
    gw=exp(-(t-(k-1)*At).^2/(4*al))/(2*sqrt(pi*al));   % Ventana de Gabor
    f_gw=f.*gw;
    G_F_GW(k,:)=abs(fft(f_gw));
end
t_new=0:At:T0;
% fd=0:1/N:1-1/N;   % Frecuencia normalizada f/fs.
fd=0:1/T0:Fs/2-1/T0;   % Frecuencia en Hz. (Medio espectro)
G_F_GW_HS=G_F_GW(:,1:N/2); % Medio espectro, hasta fs/2
H1=figure(1);
set(H1,'position',[20 30 650 650],'Menubar','none',...
        'NumberTitle','off','name', strcat('TP#5. Transformada de Gabor. Alfa = ',num2str(al)));
subplot(211)
plot(t,f)
xlabel('t(s)');
ylabel('f(t)');
grid
subplot(212)
mesh(fd,t_new,G_F_GW_HS);
axis([min(fd) max(fd) min(t_new) max(t_new) min(min(G_F_GW_HS)) max(max(G_F_GW_HS))]);
xlabel('f (Hz)');
ylabel('t (s)');
rotate3d;
H2=figure(2);
set(H2,'position',[650 20 650 650],'Menubar','none',...
        'NumberTitle','off','name', strcat('TP#5. Transformada de Gabor. Frecuencia Normalizada. Alfa = ',num2str(al)));
fdn=0:1/N:1/2-1/N;   % Frecuencia normalizada f/fs.
mesh(fdn,t_new,G_F_GW_HS);
xlabel('f/fs ()');
ylabel('t (s)');
rotate3d;