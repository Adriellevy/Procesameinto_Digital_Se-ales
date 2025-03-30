% 
% %----------------Ejercicio 1----------------
% %[t,y]=Funcion_seno_Amortiguada(10,10^-3,1,2,-0.5); 
% %[t,y]=Funcion_Sinc(10,10^-3,1,25); 
% [t,y]=Pulso(-10,10,-2,4,1);
% 
% FuncionRespuestaZ = freqz(y,1,512);
% FuncionRespuestaFourrier = fft(y);
% 
% 
% % Gráfica
% figure(1);
% subplot(2, 3,[1,4]);
% stem(t, y);
% grid on;
% xlabel('Tiempo (s)');
% ylabel('Amplitud');
% title('Señal Trigonométrica Acotada');
% 
% subplot(2, 3, 2);
% plot(abs(FuncionRespuestaZ))
% title('Módulo de la Transformada Z');
% xlabel('Re(Z)');
% ylabel('|X(z)|');
% grid on;
% subplot(2, 3, 5);
% plot((FuncionRespuestaZ))
% title('Fase de la Transformada Z');
% xlabel('Re(Z)');
% ylabel('∠X(z)');
% grid on;
% 
% subplot(2, 3, 3);
% plot(abs(FuncionRespuestaFourrier))
% title('Módulo de la Transformada Fourrier');
% xlabel('t');
% ylabel('|X(z)|');
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
% title('Señal Trigonométrica Acotada');
% 
% subplot(2,4, 2);
% plot(abs(ffta))
% title('Módulo de la T.F. sin desplazamiento');
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
% title('Señal Trigonométrica Acotada');
% 
% subplot(2,4, 4);
% plot(abs(fftb))
% title('Módulo de la T.F. con desplazamiento');
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

% Definir la secuencia na[n] en un intervalo de longitud N
N = 8; % Definir longitud de la secuencia para la DFT
n = 0:N-1;
na = [1, -2, zeros(1, N-2)]; % δ[n] - 2δ[n-1]

% Calcular la DFT de na[n]
Na = fft(na, N);

% Evitar división por cero
Na(abs(Na) < 1e-6) = 1e-6;

% Calcular la inversa en el dominio de la frecuencia
Nb = 1 ./ Na;

% Calcular la secuencia nb[n] aplicando la IDFT
nb = ifft(Nb, N);

% Verificación: Convolución de na[n] con nb[n]
conv_result = conv(na, nb, 'same');

% Graficar los resultados
figure;
subplot(3,1,1);
stem(n, real(na), 'filled');
title('Secuencia Original na[n]');
xlabel('n'); ylabel('Amplitud');

subplot(3,1,2);
stem(n, real(nb), 'filled');
title('Secuencia Inversa nb[n]');
xlabel('n'); ylabel('Amplitud');

subplot(3,1,3);
stem(n, real(conv_result), 'filled');
title('Resultado de la Convolución (Debe ser δ[n])');
xlabel('n'); ylabel('Amplitud');