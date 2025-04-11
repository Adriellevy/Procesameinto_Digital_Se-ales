% 
% %----------------Ejercicio 1----------------

% %[t,y]=Funcion_seno_Amortiguada(10,10^-3,1,2,-0.5); 
% %[t,y]=Funcion_Sinc(10,10^-3,1,25); 
% [t,y]=Pulso(0,500,0,40,1);

% %TFS freqz
% [TFS, w] = freqz(y,1,250);
% N = length(y);
% % % TFS Home made
% % N = length(y);
% % k = -N/2:N/2-1; %Divido el eje en 2 para que concuerde con tdf
% % Z = exp(1j*2*pi*k/N); % Defino la variable Z
% % TFSHomeMade = polyval(y, Z); % Evaluar la funcion en Z
% 
% %Transformada Discreta de Fourrier
% FuncionRespuestaFourrier = fft(y);
% 
% frecuencia_fft = linspace(0, 2*pi, N); % Asegurar que ambas transformadas coincidan
% 
% % Graficos
% figure(1);
% subplot(2, 3,[1,4]);
% stem(t, y);
% grid on;
% xlabel('Tiempo (s)');
% ylabel('Amplitud');
% title('Pulso');
% %axis([-80 80 -2 2]);
% subplot(2, 3, 2);
% plot(abs(TFS))
% title('Modulo de la TFS');
% xlabel('Posicion Memoria');
% ylabel('|X(z)|');
% grid on;
% 
% subplot(2, 3, 5);
% plot(angle(TFS))
% title('Fase de la TFS');
% xlabel('Posicion Memoria');
% ylabel('Fase X(z)');
% grid on;
% %axis([-4 4 -4 4]);
% 
% subplot(2, 3, 3);
% plot(abs(FuncionRespuestaFourrier(1:N/2)))
% title('Modulo de la TDF');
% xlabel('Posicion Memoria');
% ylabel('|X(w)|');
% grid on;
% subplot(2, 3, 6);
% plot(angle(FuncionRespuestaFourrier(1:N/2)))
% title('Fase de la TDF');
% xlabel('Posicion Memoria');
% ylabel('Fase X(w)');
% grid on;


%----------------Ejercicio 2----------------

% a = [1 3 5 7 5 3 1]; 

% b = [7 5 3 1 1 3 5]; 
% ffta = fft(a); 
% fftb= fft(b);
% figure(1);
% subplot(2, 4,[1,5]);
% stem(a);
% grid on;
% xlabel('Tiempo (s)');
% ylabel('Amplitud');
% title('Se帽al Trigonom茅trica Acotada');
% 
% subplot(2,4, 2);
% plot(abs(ffta))
% title('M贸dulo de la T.F. sin desplazamiento');
% xlabel('Re(Z)');
% ylabel('|X(z)|');
% grid on;
% 
% subplot(2,4, 6);
% plot(angle(ffta))
% title('Fase de la T.F. sin desplazamiento');
% xlabel('Re(Z)');
% ylabel('|X(z)|');
% grid on;
% 
% subplot(2, 4,[3,7]);
% stem(b);
% grid on;
% xlabel('Tiempo (s)');
% ylabel('Amplitud');
% title('Se帽al Trigonom茅trica Acotada');
% 
% subplot(2,4, 4);
% plot(abs(fftb))
% title('M贸dulo de la T.F. con desplazamiento');
% xlabel('Re(Z)');
% ylabel('|X(z)|');
% grid on;
% 
% subplot(2,4, 8);
% plot(angle(fftb))
% title('Fase de la T.F. con desplazamiento');
% xlabel('Re(Z)');
% ylabel('|X(z)|');
% grid on;

%----------------Ejercicio 3----------------

% % Definir la secuencia na[n] en un intervalo de longitud N
% N = 8; % Definir longitud de la secuencia para la DFT
% n = 0:N-1;
% na = [1, -2, zeros(1, N-2)]; % 未[n] - 2未[n-1]
% 
% % Calcular la DFT de na[n]
% Na = fft(na, N);
% 
% % Evitar divisi贸n por cero
% Na(abs(Na) < 1e-6) = 1e-6;
% 
% % Calcular la inversa en el dominio de la frecuencia
% Nb = 1 ./ Na;
% 
% % Calcular la secuencia nb[n] aplicando la IDFT
% nb = ifft(Nb, N);
% 
% % Verificaci贸n: Convoluci贸n de na[n] con nb[n]
% conv_result = conv(na, nb, 'same');
% 
% % Graficar los resultados
% figure;
% subplot(3,1,1);
% stem(n, real(na), 'filled');
% title('Secuencia Original na[n]');
% xlabel('n'); ylabel('Amplitud');
% 
% subplot(3,1,2);
% stem(n, real(nb), 'filled');
% title('Secuencia Inversa nb[n]');
% xlabel('n'); ylabel('Amplitud');
% 
% subplot(3,1,3);
% stem(n, real(conv_result), 'filled');
% title('Resultado de la Convoluci贸n (Debe ser 未[n])');
% xlabel('n'); ylabel('Amplitud');


%----------------Ejercicio 4----------------
%Sobre Muestreo de una funcion
% a = [2 5 3]; 
% N = length(a);
% Factor = 3;
% Func = Sobre_Muestreo(a,Factor);
% % Visualizaci髇 de resultados
% figure(1);  
% stem(0:N*Factor-1, Func, 'filled');
% legend('Sobremuestreada');
% title('Sobremuestreo de la se馻l');
%----------------Ejercicio 5----------------
% 
% a = [2 1 5 7 3 9]; 
% N = length(a);
% Factor = 2;
% Func = Sub_Muestreo(a,Factor);
% % Visualizaci髇 de resultados
% figure(2);  
% stem(1:N,a, 'b');
% hold on
% stem(1:Factor:N, Func, 'r');
% xlim([1 N]);
% title('SubMuestreo de la se馻l');

%Ejercicios 6 y 7 son fft (se ignoran)
%----------------Ejercicio 8----------------
% %A)
% A=[1 0.13 0.52 0.3] ;
% B=[0.16 -0.48 0.48 -0.16];
% [Filtroa, w]=freqz(B,A,512,'whole');
% Graficador_Freqz(Filtroa,w);
% %B)
% A=[1 0 -0.268] ;
% B=[0.634 0 -0.634];
% [Filtroa, w]=freqz(B,A,512,'whole');
% Graficador_Freqz(Filtroa,w);
% %C)
% A=[1 0 0.268] ;
% B=[0.634 0 0.634];
% [Filtroa, w]=freqz(B,A,512,'whole');
% Graficador_Freqz(Filtroa,w);
% %D)
% A=[10 -5 1] ;
% B=[1 -5 10];
% [Filtroa, w]=freqz(B,A,512,'whole');
% Graficador_Freqz(Filtroa,w);

%----------------Ejercicio 9----------------
% phi=pi/3;
% r=0.75;
% A=[1 -2*r*cos(phi) r^2] ;%yn
% B=[1];%xn
% [Filtroa, w]=freqz(B,A,512,'whole');
% Graficador_Freqz(Filtroa,w);

%----------------Ejercicio 10----------------
% %A)
% B=0.0761*[1 0 -0.7631 0 1] ;%Numerador
% A=[1 0 1.355 0 0.6196];%Denominador
% [Filtroa, w]=freqz(B,A,512,'whole');
% Graficador_Freqz(Filtroa,w);
% %B)
% B=[0.058 -0.1553 0.155 0.0518] ;%Numerador
% A=[1 1.2828 1.0388 0.341];%Denominador
% [Filtroa, w]=freqz(B,A,512,'whole');
% Graficador_Freqz(Filtroa,w);

%Ejercicio 11 y 12 realizaro por el profesor

% Leer archivo .wav
[x, fs] = audioread('introBlackDog.wav');

% Si es est閞eo, lo pasamos a mono
if size(x,2) > 1
    x = mean(x, 2);
end

% Aplicar filtrado por nivel
modo = 'amplitud';  % tambi閚 puede ser 'energia'
L = 0.01;           % Umbral (valor depende del modo)

x_filtrada = filtroPorNivel(x, modo, L);

% Graficar original y filtrada
t = (0:length(x)-1)/fs;

figure;
subplot(2,1,1);
plot(t, x);
title('Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2,1,2);
plot(t, x_filtrada);
title(['Filtrada por Nivel - Modo: ', modo]);
xlabel('Tiempo (s)');
ylabel('Amplitud');

% Escuchar se馻l original y filtrada
disp('Reproduciendo se馻l original...');
sound(x, fs);
pause(length(x)/fs + 1);

disp('Reproduciendo se馻l filtrada...');
sound(x_filtrada, fs);

