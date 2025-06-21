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
% plot(mse)

% % Automatico
N=1000;
Nro=5;
k=(0:N)';
xk=sin(2*pi*k/N);			% Señal de Entrada
dk=2*cos(2*pi*k/N);			% Señal deseada
P=1;      
w0 = zeros(1, N);       
w1 = zeros(1, N);
y = zeros(1, N);
e = zeros(1, N);  
mu = 0.5;

for i = 2:N
    y(i) = w0(i-1)*k(i) + w1(i-1)*k(i-1);    % salida del filtro
    e(i) = dk(i) - y(i);                      % error
    w0(i) = w0(i-1) + 2*mu*e(i)*k(i);        % actualizar peso 0
    w1(i) = w1(i-1) + 2*mu*e(i)*k(i-1);      % actualizar peso 1
end


[Bk,Ak,yk,ek]=ARMA_Adaptativo(xk',dk',1,0);
Bk
Ak



% Crear grilla de valores para b0 y b1
w0 = -10:0.2:10;
w1 = -10:0.2:10;
[B0, B1] = meshgrid(w0, w1);

% Inicializar matriz de error cuadrático medio
E = zeros(size(B0));

% Calcular ECM para cada par (b0, b1)
for i = 1:size(B0, 1)
    for j = 1:size(B0, 2)
        b0 = B0(i, j);
        b1 = B1(i, j);
        yk_temp = b0 * xk(2:end) + b1 * xk(1:end-1);
        ek_temp = dk(2:end) - yk_temp;
        E(i, j) = mean(ek_temp.^2);
    end
end

% Graficar superficie
figure(2);
mesh(B0, B1, E);
xlabel('b0');
ylabel('b1');
zlabel('E[ek^2]');
title('Superficie del Error Cuadrático Medio');
rotate3d on;

% Graficar curvas de nivel
figure(3);
contour(B0, B1, E, 30);
xlabel('b0');
ylabel('b1');
title('Curvas de Nivel del Error Cuadrático Medio');
grid on;
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
% [Filtro,w]= freqz([1 0 0],[1 -1.2 0.6],5000);
% figure(1)
% plot(w,abs(Filtro))
% 
% 
% N = 5000;       % Número de muestras
% mu = 0;         % Media deseada
% sigma = 0.5;      % Desvío estándar deseado
% 
% ruido = sigma * randn(N, 1)';
% ydatok=filter([1 0 0],[1 -1.2 0.6],ruido); %genero la informacion obtenida de un sensor
% 
% %Caracterizacion del sistema
% [Bk,Ak,yk,ek]=ARMA_Adaptativo(ruido',ydatok',2,2);%bk ceros, Ak polos
% [FiltroObtenido,w2]= freqz(Bk,Ak,5000); 
% Graficador_Freqz(FiltroObtenido,w2,'Caracterizacion del sistema')
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
% [Bk, Ak, yk, ek] = ARMA_Adaptativo(mic2, mic1, 15, 0);
% sound(ek,fs);
