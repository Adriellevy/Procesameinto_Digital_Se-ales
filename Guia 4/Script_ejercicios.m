%Ejercicio 1

% %Manual
% N=1000;
% k=0:0.1:100;
% frec = 10; %10hz
% w=2*pi*frec; %20pi rad/seg
% xk=sin((2*pi*k)/N);
% dk=2*cos((2*pi*k)/N);
% bk = [6 1.5]; 
% yk= bk'*xk;
% ek = yk-dk;
% mse =mean(ek.^2);
% 
%plot(mse)
% 
%Automatico
% k=0:0.1:10;
% N = length(k); % cantidad de muestras
% x = sin(2*pi*k);       % entrada
% d = 2*cos(2*pi*k); 
% w0 = zeros(1, N);       
% w1 = zeros(1, N);
% y = zeros(1, N);
% e = zeros(1, N);  
% mu = 0.5;
% 
% for k = 2:N
%     y(k) = w0(k-1)*x(k) + w1(k-1)*x(k-1);    % salida del filtro
%     e(k) = d(k) - y(k);                      % error
%     w0(k) = w0(k-1) + 2*mu*e(k)*x(k);        % actualizar peso 0
%     w1(k) = w1(k-1) + 2*mu*e(k)*x(k-1);      % actualizar peso 1
% end
% 
% figure(1)
% subplot(3,1,1)
% plot(w0); hold on; plot(w1); legend('w0','w1');
% title('Evolucion de los pesos');
% hold off;
% subplot(3,1,2)
% plot(d); hold on; plot(y); legend('d','y');
% title('Comparacion funciones');
% hold off;
% subplot(3,1,3)
% contour(w0,w1,ek,30)
% title('Comparacion funciones');
%Ver como encarar el tema de conseguir los autovectores

%Ejercicio 2 --> hacer la identificacion del sistema
% N=512;
% fp=60; %hz 
% wc = 2*pi*fp;
% fs=1000; %hz
% Ts=1/fs;
% nTs = 0:Ts:(N-1)*Ts;
% xk = 2*sin(2*pi*50*nTs)+sin(2*pi*300*nTs); 
% dk = 2*sin(2*pi*50*nTs);
% 
% [numz,denz]=bilinear(wc,[1 wc],fs);
% ydk=filter(numz,denz,xk);
% [Bk,Ak,yk,ek]=ARMA_Adaptativo(xk',ydk',1,1);
% subplot(4,1,1) 
% plot(xk) 
% subplot(4,1,2)
% plot(ydk)
% subplot(4,1,3) 
% plot(yk)
% subplot(4,1,4) 
% plot(ek)
% [H,w]=freqz(Bk,Ak,512);
% Graficador_Freqz(H,w,'Grafico');
% hold off;
% [num_s, den_s] = invbilinear(Bk',Ak', Ts);
% figure(3)
% hs = tf(num_s,den_s);
% bode(hs) %La frecuencia de corte da 377

%Ejercicio 3
% Tengo que eliminar el ruido (hay que meterle ruido a una señal
%estocastica) 

% 
% N = 500;       % Número de muestras
% mu = 5;         % Media deseada
% sigma = 5;      % Desvío estándar deseado
% 
% ruido = sigma * randn(N, 1)';
% fp=60; %hz 
% wc = 2*pi*fp;
% fs=1000; %hz
% Ts=1/fs;
% nTs = 0:Ts:(N-1)*Ts;
% xk = 2*sin(2*pi*50*nTs)+ruido; 
% dk = 2*sin(2*pi*50*nTs);
% [numz,denz]=bilinear(wc,[1 wc],fs);
% ydk=filter(numz,denz,xk);
% [Bk,Ak,yk,ek]=ARMA_Adaptativo(xk',ydk',1,1);
% subplot(4,1,1) 
% plot(xk) 
% subplot(4,1,2)
% plot(ydk)
% subplot(4,1,3) 
% plot(yk)
% subplot(4,1,4) 
% plot(ek)
% [H,w]=freqz(Bk,Ak,512);


%Ejercicio 4 en este ejerccio se entiende que se puede obtener la
%caracterizacion del sistema (esto es el diagnostico)
% 
% N = 500;       % Número de muestras
% mu = 5;         % Media deseada
% sigma = 5;      % Desvío estándar deseado
% 
% ruido = sigma * randn(N, 1)';
% ydatok=filter([1 0 0],[1 -1.2 0.6],ruido); %genero la informacion obtenida de un sensor
% [Bk,Ak,yk,ek]=ARMA_Adaptativo(xk',ydk',0,2);
% Bk 
% Ak

%Ejercicio 5 
% % Abrir el archivo del ecg
% fid = fopen('ecg_nt_Ejercicio_5.txt', 'r');
% if fid == -1
%     error('No se pudo abrir el archivo.');
% end
% 
% % Leer datos (asume una columna de datos numéricos)
% ecg = fscanf(fid, '%f\n');
% 
% % Cerrar archivo
% fclose(fid);
% N=550;
% fs=250; %hz
% Ts=1/fs;
% nTs = 0:Ts:(N-1)*Ts;
% ruido = sin(2*pi*50*nTs);
% [Bk,Ak,yk,ek]=ARMA_Adaptativo(ruido',ecg',1,1);
% 
% subplot(4,1,1) 
% plot(ecg) 
% subplot(4,1,2)
% plot(ruido)
% subplot(4,1,3) 
% plot(yk)
% subplot(4,1,4) 
% plot(ek)
%Ejercicio 6 predictivo
%voy a tener que calcular el orden del filtro y 
% fid = fopen('ecg_nt_long_E_6.txt', 'r');
% if fid == -1
%     error('No se pudo abrir el archivo.');
% end
% 
% % Leer datos (asume una columna de datos numéricos)
% ecg = fscanf(fid, '%f\n');
% 
% % Cerrar archivo
% fclose(fid);
% shift = 50;
% ecgdesplazado = [zeros(shift,1); ecg(1:end-shift)]; %desplazamiento de
% tiempo me ayuda a predecir el electrocardiograma
% [Bk,Ak,yk,ek]=ARMA_Adaptativo(ecgdesplazado,ecg,15,0);
% 
% subplot(4,1,1) 
% plot(ecg) 
% subplot(4,1,2)
% plot(ecgdesplazado)
% subplot(4,1,3) 
% plot(yk)
% subplot(4,1,4) 
% plot(ek(500:1000))


%Ejercicio 7 

%Ejercicio 8

% [Bk,Ak,yk,ek]=ARMA_Adaptativo(mic2,mic1,15,0);
% plot(ek)
%[Bk, Ak, yk, ek] = ARMA_Adaptativo(mic2, mic1, 15, 0);
sound(ek(250000:280000),fs);
