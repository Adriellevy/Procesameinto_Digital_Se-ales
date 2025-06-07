%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Cálculo de la Transformada Ondita Discreta Inversa
%%%
%%% function Mat_s=IDWT_FMP_Second_Way(Mat,L,lp)
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
%%%   Mat_s: Matriz de detalles y Aproximación 
%%%         Fila1:       Det 1
%%%         Fila2:       Det 2
%%%            :             :
%%%         FilaNivel:   Det Nivel
%%%         FilaNivel+1: App Nivel 
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

function Mat_s=IDWT_FMP_Second_Way_Ejer_14(Mat,L,lp)

if nargin<3
    Return;
end
N=length(lp); % Longitud de los filtros de descomposición y análisis
n=0:N-1;
l=lp(N-n);  % Filtro Pasa Bajos de descomposición
hp=(-1).^n.*l; % Filtro Pasa Altos de reconstrucción
Niv=length(L)-1;
N_sig=sum(L);  % Longitud de la señal desccompuesta
Mat_s=zeros(Niv+1,N_sig); % Matriz que contiene los detalles y las aproximación final
%Det_i=Mat(Niv+1,1:L(Niv+1));
k=1;
Cd_i=Mat(k,1:L(1));
Cdi_up=dyadup(Cd_i);
Cdi_up=Cdi_up(2:length(Cdi_up)-1); % Quitado de ceros al inicio y fin ya que dyadup los agrega
Det_1=conv(Cdi_up,hp/sqrt(2));  % Detalle 1
Mat_s(k,:)=Det_1;  % Guardado del Detalle 1
if Niv>1
    for k=2:Niv
        Cd_i=Mat(k,1:L(k));  % Coeficientes del detalle i-ésimo
        Cdi_up=dyadup(Cd_i); % Sobremuestreo del detalle i-ésimo
        Cdi_up=Cdi_up(2:length(Cdi_up)-1); % Quitado de ceros al inicio y fin ya que dyadup los agrega
        Det_i=conv(Cdi_up,hp/sqrt(2));
        for j=1:k-1
            Det_i=dyadup(Det_i);
            Det_i=Det_i(2:length(Det_i)-1);
            Det_i=conv(Det_i,lp/sqrt(2));
        end
        Mat_s(k,:)=Det_i;  % Guardado del detalle i-ésimo
        if k==Niv
            Ca_i=Mat(k+1,1:L(k+1));  % Coeficientes de la última aproximación
            Cai_up=dyadup(Ca_i);   % Sobremuestreo de la aproximación
            Cai_up=Cai_up(2:length(Cai_up)-1); % Quitado de ceros al inicio y fin ya que dyadup los agrega
            App_i=conv(Cai_up,lp/sqrt(2));
            for j=1:k-1
                App_i=dyadup(App_i);
                App_i=App_i(2:length(App_i)-1);
                App_i=conv(App_i,lp/sqrt(2));
            end
            Mat_s(k+1,:)=App_i;    
        end
    end
else
    Ca_i=Mat(k+1,1:L(k+1));  % Coeficientes de la última aproximación
    Cai_up=dyadup(Ca_i);   % Sobremuestreo de la aproximación
    Cai_up=Cai_up(2:length(Cai_up)-1); % Quitado de ceros al inicio y fin ya que dyadup los agrega
    App_i=conv(Cai_up,lp/sqrt(2));
    for j=1:k-1
        App_i=dyadup(App_i);
        App_i=App_i(2:length(App_i)-1);
        App_i=conv(App_i,lp/sqrt(2));
    end
    Mat_s(Niv+1,:)=App_i;    
 end