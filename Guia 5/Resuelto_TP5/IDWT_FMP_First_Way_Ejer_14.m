%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Cálculo de la Transformada Ondita Discreta Inversa
%%%
%%% function s=IDWT_FMP_First_Way(Mat,L,lp)
%%%
%%% Parámetros de Entrada:
%%%
%%%    Mat: Matriz con el siguiente contenido:
%%%         Fila1:       Coef Det 1
%%%         Fila2:       Coef Det 2
%%%            :             :
%%%         FilaNivel:   Coef Det Nivel
%%%         FilaNivel+1: Coef App Nivel 
%%%    L: Vector conteniendo la longitud de los coef de det y app
%%%         L(1):       Long Coef Det 1
%%%         L(2):       Long Coef Det 2
%%%           :                :
%%%         L(Nivel):   Long Coef Det Nivel
%%%         L(Nivel+1): Long Coef App Nivel
%%%    lp: Filtro Pasa Bajos de Reconstrucción
%%%
%%% Parámetros de Salida
%%%
%%%   s: Señal reconstruida 
%%%
%%% IMPORTANTE:
%%% hp, filtro pasa bajos de reconstrucción, no tiene que estar
%%% normalizado, ya que dentro de esta rutina se realiza el escalamiento de
%%% los filtros por sqrt(2) para el escalamiento de energías
%%%
%%% Versión Propia del comando MatLab waverec 
%%%
%%% Procesamiento Digital de Señales
%%%
%%% Dr. Ing. Franco Martin Pessana
%%% 
%%% Facultad de Ingeniería y Ciencias Exactas y Naturales
%%% Universidad Favaloro
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function s=IDWT_FMP_First_Way_Ejer_14(Mat,L,lp)

if nargin<3
    Return;
end

N=length(lp); % Longitud de los filtros de descomposición y análisis
n=0:N-1;
l=lp(N-n);  % Filtro Pasa Bajos de descomposición
hp=(-1).^n.*l; % Filtro Pasa Altos de reconstrucción
Niv=length(L)-1;
N_sig=sum(L);
Ca_i=Mat(Niv+1,1:L(Niv+1));
for k=1:Niv
    Cd_i=Mat(Niv+1-k,1:L(Niv+1-k));
    Cdi_up=dyadup(Cd_i);
    Cdi_up=Cdi_up(2:length(Cdi_up)-1);
    Cai_up=dyadup(Ca_i);
    Cai_up=Cai_up(2:length(Cai_up)-1);
    Ca_i=conv(Cdi_up,hp/sqrt(2))+conv(Cai_up,lp/sqrt(2));
end
s=wkeep(Ca_i,N_sig);