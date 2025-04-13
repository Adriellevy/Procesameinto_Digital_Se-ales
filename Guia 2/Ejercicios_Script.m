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
% title('SeÃ±al TrigonomÃ©trica Acotada');
% 
% subplot(2,4, 2);
% plot(abs(ffta))
% title('MÃ³dulo de la T.F. sin desplazamiento');
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
% title('SeÃ±al TrigonomÃ©trica Acotada');
% 
% subplot(2,4, 4);
% plot(abs(fftb))
% title('MÃ³dulo de la T.F. con desplazamiento');
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
% na = [1, -2, zeros(1, N-2)]; % Î´[n] - 2Î´[n-1]
% 
% % Calcular la DFT de na[n]
% Na = fft(na, N);
% 
% % Evitar divisiÃ³n por cero
% Na(abs(Na) < 1e-6) = 1e-6;
% 
% % Calcular la inversa en el dominio de la frecuencia
% Nb = 1 ./ Na;
% 
% % Calcular la secuencia nb[n] aplicando la IDFT
% nb = ifft(Nb, N);
% 
% % VerificaciÃ³n: ConvoluciÃ³n de na[n] con nb[n]
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
% title('Resultado de la ConvoluciÃ³n (Debe ser Î´[n])');
% xlabel('n'); ylabel('Amplitud');


%----------------Ejercicio 4----------------
%Sobre Muestreo de una funcion
% a = [2 5 3]; 
% N = length(a);
% Factor = 3;
% Func = Sobre_Muestreo(a,Factor);
% % Visualización de resultados
% figure(1);  
% stem(0:N*Factor-1, Func, 'filled');
% legend('Sobremuestreada');
% title('Sobremuestreo de la señal');
%----------------Ejercicio 5----------------
% 
% a = [2 1 5 7 3 9]; 
% N = length(a);
% Factor = 2;
% Func = Sub_Muestreo(a,Factor);
% % Visualización de resultados
% figure(2);  
% stem(1:N,a, 'b');
% hold on
% stem(1:Factor:N, Func, 'r');
% xlim([1 N]);
% title('SubMuestreo de la señal');

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

%Ejercicio 11 realizado por el profesor

%Ejercicio 12 y 13 Retardo de Grupo, no me interesa

%Ejercicio 14 y 15 Grafico de filtros y funciones

%Ejercicio 16 y 17 
% Leer archivo .wav
% [x, fs] = audioread('introBlackDog.wav');
% 
% % Si es estéreo, lo pasamos a mono
% if size(x,2) > 1
%     x = mean(x, 2);
% end
% Modo='energia';
% Valor=50;
% x_filtrada = FiltroPorNivel(x,Valor, Modo);
% 
% % Graficar original y filtrada
% t = (0:length(x)-1)/fs;
% 
% figure(3);
% subplot(2,1,1);
% plot(t, x);
% title('Original');
% xlabel('Tiempo (s)');
% ylabel('Amplitud');
% 
% subplot(2,1,2);
% plot(t, x_filtrada);
% title(sprintf('Filtrada por %s y valor: %0.f', Modo, Valor))
% xlabel('Tiempo (s)');
% ylabel('Amplitud');
% 
% % 
% % title('Fourrier');
% % xlabel('Memoria');
% % ylabel('|W(s)|');
% % Escuchar señal original y filtrada
% % disp('Reproduciendo señal original...');
% % sound(x, fs);
% % pause(length(x)/fs + 1);
% 
% disp('Reproduciendo señal filtrada...');
% sound(x_filtrada, fs);

%_________________APLICACIONES________________________

%----------------Ejercicio 18 submuestreo---------------- 

% [t,x]=Pulso(0,100,0,50,1,10^-1);
% XSub_muestreada = Sub_Muestreo(x,100); 
% 
% figure(8)
% stem(t,x)
% 
% figure(9)
% stem(t(1:100:end),XSub_muestreada)
% [t,x]=Pulso(0,100,0,50,1,10^-3);
% TDF = fft(x);
% Graficador_Fft(TDF)

%----------------Ejercicio 19 submuestreo----------------
% [t,x] = Pulso(0, 100, 0, 50, 1, 10^-1);
% M_vals = [2, 4, 8];
% figure(3);
% for i = 1:length(M_vals)

%     M = M_vals(i);
% 
%     % Submuestreo
%     xM = Sub_Muestreo(x, M);
%     tM = t(1:M:end);
% 
%     % Espectro de Fourier
%     N = length(xM);
%     f = (-N/2:N/2-1)*(1/(N*(M*10^-3))); % Escala de frecuencias correcta para Ts*M
%     X_f = fftshift(abs(fft(xM)));
% 
%     subplot(length(M_vals),2,2*i-1);
%     stem(tM, xM, 'filled');
%     title(['Señal Submuestreada M = ', num2str(M)]);
%     xlabel('Tiempo [s]');
%     ylabel('Amplitud');
% 
%     subplot(length(M_vals),2,2*i);
%     plot(f, X_f);
%     title(['Espectro de Fourier M = ', num2str(M)]);
%     xlabel('Frecuencia [Hz]');
%     ylabel('|X(f)|');
% 
% end

% ----------------Ejercicio 20 submuestreo----------------
% [t,x_original]=Pulso(0,100,0,50,1,10^-1);
% N = length(x_original);
% f = (-N/2:N/2-1)*(1/(N*10^-3));
% X = fftshift(fft(x_original));
% posicionesSubplot= [1 3 5];
% figure(4);
% for i = 1:length(M_vals)
%     M = M_vals(i);
%     wc = round(N / (2*M));  % índice correspondiente a w_c = pi/M
% 
%     % Crear filtro en frecuencia con tu función
%     H = FiltroPasaBajos(N, wc);  % H está centrado en 0 por definición
% 
%     % Aplicar filtro sobre el espectro
%     X_filtered = X .* H;
% 
%     % Volver al dominio temporal
%     x_filtered = fftshift(ifft(X_filtered));
% 
%     subplot(3,2,posicionesSubplot(i))
%     plot(t, H);
%     subplot(3,2,2*i)
%     plot(t, x_filtered);
%     title(['Señal Filtrada Idealmente para M = ', num2str(M), ' (w_c = \pi/', num2str(M), ')']);
%     xlabel('Tiempo [s]');
%     ylabel('Amplitud');
% end

%-----------------Ejercicio 21 submuestreo----------------
%Quiere comparar ejercicio 18 con 19 y 20


%-----------------Ejercicio 22 y 23 sobremuestreo----------------
% Señal original
% [t, x] = Funcion_Sinc(5, 5e-2, 1, 0); 
% 
% figure('Name','Ejercicio 22 y 23')
% subplot(3,2,1)
% stem(t, x, 'filled');
% title('Señal Original: x[n]');
% xlabel('Tiempo [s]');
% ylabel('Amplitud');
% grid on;
% xlim([min(t) max(t)]);
% 
% % Aumento de resolución
% Factor = 3; 
% x_aumentada = Sobre_Muestreo(x, Factor);
% 
% % Transformada Discreta de Fourier (TDF)
% TDF = fft(x_aumentada);
% 
% M = length(TDF);
% w = linspace(0, 2*pi, M); % Crear eje de frecuencia normalizada
% 
% % Señal sobremuestreada
% subplot(3,2,2)
% stem(x_aumentada, 'filled');
% title(['Señal Sobremuestreada (Factor = ', num2str(Factor), ')']);
% xlabel('n');
% ylabel('Amplitud');
% grid on;
% xlim([0 length(x_aumentada)]);
% 
% 
% % Espectro de la señal sobremuestreada
% subplot(3,2,3)
% plot(w, abs(TDF), 'LineWidth', 1.5);
% title('Magnitud de la FFT de x[n]');
% xlabel('Frecuencia Normalizada \omega [rad/muestra]');
% ylabel('|X(e^{j\omega})|');
% grid on;
% xlim([0 2*pi]);
% 
% % Filtro pasa bajos
% N = length(TDF);
% wc = round(N / (2 * Factor)); 
% Filtro = FiltroPasaBajos(N, wc);
% 
% subplot(3,2,4)
% stem(w, Filtro, 'LineWidth', 1.5);
% title(['Filtro Pasa Bajos (wc = ', num2str(wc), ')']);
% xlabel('Frecuencia Normalizada \omega [rad/muestra]');
% ylabel('H(e^{j\omega})');
% grid on;
% xlim([0 2*pi]);
% 
% % Aplicación del filtro en frecuencia
% Resultado = TDF .* Filtro;
% 
% subplot(3,2,5)
% plot(w, abs(Resultado), 'LineWidth', 1.5);
% title(['FFT Filtrada: |X[k]·H[k]|']);
% xlabel('Frecuencia Normalizada \omega [rad/muestra]');
% ylabel('|Y(e^{j\omega})|');
% grid on;
% xlim([0 2*pi]);
% 
% % Señal filtrada en el dominio temporal
% v_filtrada = ifft(Resultado);
% t_filtrada = linspace(0,5,length(v_filtrada));
% 
% subplot(3,2,6)
% stem(t_filtrada,real(v_filtrada), 'filled');
% title('Señal Filtrada en el Tiempo');
% xlabel('n');
% ylabel('Amplitud');
% grid on;


%-----------------Ejercicio 24 correlacion---------------
% 
% [t,x]=Pulso(0,100,0,10,1,10^-1);
% rx = Correlacion(x);        % Autocorrelación de x
% 
% figure(3)
% subplot(2,1,1); 
% stem(x); 
% title('Pulso');
% 
% subplot(2,1,2);  
% t_nueva= 1:length(rx)/2;
% x_obtenida=rx(1:length(rx)/2)/max(abs(rx));
% stem(t_nueva,x_obtenida); 
% title('Autocorrelacion Pulso');
% 

%-----------------Ejercicio 25 sonar---------------


%-----------------Ejercicio 25 sonar---------------

%-----------------Ejercicio 26 Ventanas de visualizacion---------------
%-----------------Ejercicio 27 Ventana cuadrada------------------------
% Fs = 10000;         % Frecuencia de muestreo
% N = 1024;           % Tamaño de la FFT
% A1 = 1;
% A2 = 0.75;
% 
% L_vals = [64, 512]; % Longitudes de ventana
% w1_list = [2*pi/6,  2*pi/14, 2*pi/14, 2*pi/14];
% w2_list = [2*pi/3,  4*pi/15, 2*pi/12, 2*pi/12];
% L_list  = [64,      64,      64,      512];
% 
% for k = 1:4
%     w1 = w1_list(k);
%     w2 = w2_list(k);
%     L = L_list(k);
% 
%     n = 0:L-1;
%     x = A1*cos(w1*n) + A2*cos(w2*n);
% 
%     % Aplicar ventana rectangular
%     ventana = ones(1, L);
%     x_win = x .* ventana;
% 
%     % Aplicar filtro pasa bajos (usando tu función)
%     wc = floor(L/8); % ejemplo: eliminar centro
%     Filtro = FiltroPasaBajos(L, wc);
%     x_filt = x_win .* Filtro;
% 
%     % FFT
%     X = fft(x_filt, N);
%     f = linspace(0, Fs, N);
% 
%     % Plot
%     figure;
%     plot(f, abs(X));
%     title(sprintf('Espectro - Caso %d: L = %d', k, L));
%     xlabel('Frecuencia (Hz)');
%     ylabel('|X(f)|');
%     grid on;
% end


