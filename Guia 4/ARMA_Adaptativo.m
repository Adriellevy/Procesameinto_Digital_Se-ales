%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Rutina que permite el calculo de filtros adaptativos IIR (o FIR)
%%% mediante la sucesion de 2 filtros FIR. Utiliza el metodo del
%%% descenso mas escalonado junto con LMS
%%% 
%%%  Sintaxis:
%%%
%%%    [B,A,yk,ek]=Arma_Adaptativo(x,d,P,Q)
%%%          
%%%  Donde:
%%%
%%% Parámetros de Salida:
%%%          B:  Coeficientes del numerador del Filtro Adaptativo IIR
%%%          A:  Coeficientes del denominador del Filtro Adaptativo IIR
%%%          yk: salida del filtro adaptativo IIR
%%%          ek: Señal de error en el proceso de adapatacion
%%%
%%% Parámetros de entrada:
%%%          x: Señal a la entrada del filtro digital adaptativo (Columna)
%%%          d: Señal deseada del proceso de adapatacion (Columna)
%%%          P: Orden del Numerador (Cantidad de ceros)
%%%          Q: Orden del Denominador (Cantidad de Polos)
%%%
%%% Dr. Ing. Franco Martin Pessana
%%% Facultad de Ingeniería y Ciencias Exactas y Naturales
%%% Universidad Favaloro
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B,A,yk,ek]=ARMA_Adaptativo(x,d,P,Q)

N=length(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Sistema Digital de ARMA(P,Q) o AR(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pot_x=mean(x.*x);
Pot_d=mean(d.*d);
mu1=1/(20*(P+1)*Pot_x);	% vel. de Convergencia del Algoritmo para ceros
mu2=1/(20*(Q+1)*Pot_d);	% vel. de Convergencia del Algoritmo para polos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Recordar: Cuanto más pequeña la señal a adaptar, mas grande debe ser mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yk=zeros(size(x));  % Señal de salida inicial del Filtro IIR
ek=zeros(size(x));	% Señal inicial de error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=max([P Q]);
if Q==0 % Se está trabajando con un FIR Adaptativo
    Bk=zeros(P+1,1);	% Vector Inicial de Pesos del FIR de los ceros
    for i=P+1:N
        x_i=x(i:-1:i-P);  % Vector de entrada x i-ésimo
        yk(i)=Bk'*x_i;	        % señal a la salida del FIR
        ek(i)=d(i)-yk(i);		    % Señal de error
        Bk=Bk+2*mu1*ek(i)*x_i;	    % Vector de pesos i-ésimo de los ceros
    end
    B=Bk;
    A=1;
else % Se está trabajando con un IIR Adaptativo
    Bk=zeros(P+1,1);	% Vector Inicial de Psos del FIR de los ceros
    Ak=zeros(Q,1);	% Vector Inicial de Pesos del FIR de los polos
    for i=M+1:N-1
        x_i=x(i:-1:i-P);  % Vector x i-ésimo
        d_i=d(i-1:-1:i-Q);
        yk(i)=Bk'*x_i+Ak'*d_i;	% señal a la salida del IIR
        ek(i)=d(i)-yk(i);		    % Señal de error
        Bk=Bk+2*mu1*ek(i)*x_i;	% Vector de pesos i-ésimo de los ceros
        Ak=Ak+2*mu2*ek(i)*d_i;	% Vector de pesos i-ésimo de los polos
    end
    B=Bk;
    A=[1;-Ak];
end