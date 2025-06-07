%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%% Ejercicio Nro. 3)
%%%
%%% Unidad N° 5. Transformada Ondita
%%%
%%% Propiedad importante de sobremuestreo y Ts=1
%%%
%%% Procesamiento Digital de Señales
%%%
%%% Dr. Ing. Franco Martin Pessana
%%% 
%%% Facultad de Ingeniería y Ciencias Exactas y Naturales
%%% Universidad Favaloro
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Ejercicio_3(a)

if nargin<1
    a=4;
end
t1=-5*a:5*a;
tdesp=0;
al=1/2;
gw1=(exp(-(t1-tdesp).^2/(4*al*a^2))/(2*sqrt(pi*al)));   % Ventana de Gabor. Primera forma
t2=-5:1/a:5;
gw2=(exp(-(t2-tdesp).^2/(4*al))/(2*sqrt(pi*al)));   % Ventana de Gabor. Segunda forma
H1=figure(1);
set(H1,'position',[20 30 1200 700],'Menubar','none',...
        'NumberTitle','off','name','Invarianza al escalamiento de w(t/a) con intervalo [-5a,5a] y g(t) con intervalo [-5,5]');
subplot(221)
stem(gw1);
xlabel('n ()');
ylabel('w(t/a) con intervalo [-5a,5a]');
grid
subplot(223)
stem(gw2);
xlabel('n ()');
ylabel('w(t) con intervalo [-5,5]');
grid
subplot(222)
stem(t1,gw1);
xlabel('t (s)');
ylabel('w(t/a) con intervalo [-5a,5a]');
grid
subplot(224)
stem(t2,gw2);
xlabel('t (s)');
ylabel('w(t) con intervalo [-5,5]');
grid