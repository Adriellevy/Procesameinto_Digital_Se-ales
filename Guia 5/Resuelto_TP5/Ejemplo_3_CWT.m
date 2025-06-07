%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Ejemplo de Cálculo de CWT chirp signal (senoide con incremento de
%%% frecuencia, utilizando Ondita Morlet
%%%
%%% MSc Ing. Franco Martin Pessana
%%% Procesamiento Avanzado de Señales
%%% Facultad Regional Buenos Aires
%%% Universidad Tecnológica Nacional
%%% Doctorado en Ingeniería
%%% Ciclo Lectivo 2012
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X=Ejemplo_3_CWT

fs=1000;
Ts=1/fs;
Tf=2;
Fi=0;
Ff=150;
t=0:Ts:Tf-Ts;                    
x=chirp(t,Fi,Tf,Ff);
H1=figure(1);
set(H1,'NumberTitle','off','Menubar','none','name','Ejemplo #3. Ondita Continua. Chirp Signal','position',[60 30 1000 700]);
subplot(121)
plot(t,x);
xlabel('t');
ylabel('x(t)');
title('Chirp signal')
grid;
subplot(122)
Esc=1:64;
X=cwt(x,Esc,'morl','plot');
Y=abs(round(X)); % Vista de la matriz como en el libro para graficar
H2=figure(2);
set(H2,'NumberTitle','off','Menubar','none','name','Ejemplo #3. Vista 3D','position',[100 100 600 600]);
mesh(abs(X));
xlabel('n');
ylabel('esc');
title('|X(b,a)|')
rotate3d;