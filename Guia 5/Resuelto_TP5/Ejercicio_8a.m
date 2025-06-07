%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%% Ejercicio Nro. 8a) de la Guía_5
%%%
%%% Unidad N° 5
%%%
%%% Armando del Banco de Filtros de Aproximación
%%% 
%%% funcion padre fi_sinc(t)
%%% 
%%% IMPORTANTE:
%%% El eje de frecuencias normalizado f/fs va hasta 1, pero aquí esta
%%% representadala mitad del espectro,
%%% es por ello, que el primer filtro que corta en pi/2, equivalente a fs/4
%%% y normalizado en 0.25, está indicado aquí como 0.5 
%%%
%%% Procesamiento Digital de Señales
%%%
%%% Dr. Ing. Franco Martin Pessana
%%% 
%%% Facultad de Ingeniería y Ciencias Exactas y Naturales
%%% Universidad Favaloro
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ejercicio_8a(M)

if nargin<1
    M=6;
end
t=-50*2^M:50*2^M; % escala de tiempos adimensional
N=length(t);
fn=(0:1/N:1-1/N); % escala normalizada de frecuencia
m_fi=zeros(M,N); % Matriz que contiene las onditas padre escaladas
M_FI=zeros(M,N); % Matriz que contiene el espectro de las onditas escaladas
for j=1:M
    m_fi(j,:)=2^(-j+1)*fi_sinc(2^(-j+1)*t); %%%%% feval('fi_sinc',2^(-j+1)*t);
    M_FI(j,:)=abs(fft(m_fi(j,:)));
end
H1=figure(1);
set(H1,'position',[20 30 1200 700],'Menubar','none',...
        'NumberTitle','off','name',strcat('TP5. Escalamiento Ondita Padre. Tiempo. Cantidad= ',num2str(M), '. Graficado de una sola con j=',num2str(M)));
plot(t,m_fi(M,:));
xlabel('t ()');
ylabel(strcat('Ondita Padre Fi(t) con j = ',num2str(M)));
grid
H2=figure(2);
set(H2,'position',[60 30 1200 700],'Menubar','none',...
        'NumberTitle','off','name',strcat('TP5. Escalamiento Ondita Padre. Frecuencia. Cantidad= ',num2str(M)));
%plot(fn,M_W);
plot(fn(1:round(N/2)),M_FI(:,1:round(N/2)));
xlabel('f/fs ()');
ylabel('H(f) escaladas y normalizadas');
grid;