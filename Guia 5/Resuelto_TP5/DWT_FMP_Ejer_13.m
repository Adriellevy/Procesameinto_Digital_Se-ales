%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Cálculo de la Transformada Ondita Discreta
%%%
%%% function [Mat,L,Coef]=DWT_FMP_Ejer_13(s,lp,Niv)
%%%
%%% Parámetros de Entrada:
%%%
%%%    s:   Señal a Transformar
%%%    lp:  FPBj de Reconstrucción
%%%    Niv: Nivel de descomposición
%%%
%%% Parámetros de Salida:
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
%%%    Coef: Vector con concatenación de Coef (como MatLab)
%%%    Coef=[Coef_Det_1|Coef_Det_2| ... |Coef_Det_Nivel||Coef_App_Nivel]
%%%
%%% IMPORTANTE:
%%% hp, filtro pasa bajos de reconstrucción, no tiene que estar
%%% normalizado, ya que dentro de esta rutina se realiza el escalamiento de
%%% los filtros por sqrt(2) para el escalamiento de energías
%%%
%%% Versión Propia del comando MatLab wavedec
%%%
%%% Procesamiento Digital de Señales
%%%
%%% Dr. Ing. Franco Martin Pessana
%%% 
%%% Facultad de Ingeniería y Ciencias Exactas y Naturales
%%% Universidad Favaloro
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Mat,L,Coef]=DWT_FMP_Ejer_13(s,lp,Niv)

if nargin<3
    Niv=1;
end
N=length(lp); % Longitud de los filtros de descomposición y análisis
V=floor(log2(length(s)));
if Niv>V
    disp('No puede haber mas niveles que log2(longitud señal). Ingrese menos niveles');
    return;
end
n=0:N-1;
l=lp(N-n);  % Filtro Pasa Bajos de descomposición
hp=(-1).^n.*l; % Filtro Pasa Altos de reconstrucción
h=hp(N-n); % Filtro Pasa Altos de Descomposición
Mat=zeros(Niv+1,floor((length(s)+length(h)-1)/2)); % Matriz que contiene coeficientes de detalles y aproximación
L=zeros(1,Niv+1); % Vector que contiene las longitudes de los coeficientes de detalle y aproximación
Coef=zeros(1,length(s)); % Vector que contiene los coeficientes de detalle y aproximación
CA_i=s; % Coeficiente de aproximación de Nivel 0 (la señal de entrada)
for k=1:Niv
    CD_i=dyaddown(conv(CA_i,h/sqrt(2)));  % Coeficiente de detalle i-ésimo. El submuestreo deja los componentes impares
    CA_i=dyaddown(conv(CA_i,l/sqrt(2)));  % Coeficiente de aproximación i-ésimo. Con downsample se dejan las impares
    L(k)=length(CD_i);
    Coef(1+sum(L(1:k-1)):sum(L(1:k)))=CD_i;
    Mat(k,1:length(CD_i))=CD_i;
    if k==Niv
        L(k+1)=length(CA_i);
        Coef(1+sum(L(1:k)):sum(L(1:k+1)))=CA_i;
        Mat(k+1,1:length(CA_i))=CA_i;
    end    
end
% [C,Long]=wavedec(s,Niv,'haar');