%-------------------Ejercicio 1---------------------
% a = [1, -0.8741, 0.9217, -0.2672];
% b = [0.1866, 0.2036, 0.2036, 0.1866];
% 
% % Representaci贸n en espacio de estados
% [A, B, C, D] = tf2ss(b, a);
% 
% % Expansi贸n en fracciones simples
% [r, p, k] = residuez(b, a);
% 
% [H, w] = freqz(b, a, 512);
% figure;
% subplot(2,1,1);
% plot(w, abs(H));
% title('Respuesta en Frecuencia - Magnitud');
% xlabel('\omega');
% ylabel('|H(e^{j\omega})|');
% 
% subplot(2,1,2);
% plot(w, angle(H));
% title('Respuesta en Frecuencia - Fase');
% xlabel('\omega');
% ylabel('Fase (rad)');
% 
% imp = [1; zeros(99,1)];
% h = filter(b, a, imp);
% 
% % Comparaci贸n con fracciones simples
% h_analitica = zeros(100, 1);
% for i = 1:length(r)
%     h_analitica = h_analitica + real(r(i) * (p(i).^(0:99)'));
% end
% 
% % Comparaci贸n visual
% figure;
% plot(0:99, h, 'b', 0:99, h_analitica, 'r--');
% legend('Respuesta filter', 'Respuesta fracciones simples');
% title('Comparaci贸n respuesta al impulso');

%-------------------Ejercicio 2---------------------
% Datos del problema
% p = [0.9, 0.6718 + 0.6718j, 0.6718 - 0.6718j];
% z = [1, 1j, -1j];
% k = 0.771;
% 
% [b, a] = zp2tf_custom(z, p, k);
% 
% v1 = ones(100,1);
% v2 = ones(100,1) * -1; 
% v2(1:6) = 1; % [1 1 1 1 1 1 -1 -1 ...]
% v3 = repmat([1; 0], 50, 1);         % [1 0 1 0 ... 100 elementos]
% 
% y1 = filter(b, a, v1);
% y2 = filter(b, a, v2);
% y3 = filter(b, a, v3);
% 
% v = repmat(v1, 25, 1);  % donde vp es una secuencia de longitud 4
% y = filter(b, a, v);
%-------------------Ejercicio 3---------------------
% Datos del problema
% A = [0     0.8580     0;
%      1.3639 0.5111 0.8580;
%      0.3639     0     0];
% b = [0; 1; 1];
% c = [1.3629 0.6019 0.3074];
% d = 1;
% 
% [b_tf, a_tf] = ss2tf(A, b, c, d);
% polos = roots(a_tf);
% autovalores = eig(A);  % deber铆an coincidir
% N = 50;
% x = zeros(3, N+1);
% h = zeros(1, N+1);
% v = [1, zeros(1,N)];
% 
% for n = 1:N+1
%     h(n) = c * x(:,n) + d * v(n);
%     if n < N+1
%         x(:,n+1) = A * x(:,n) + b * v(n);
%     end
% end
% 
% figure;
% subplot(2,2,1); plot(0:N, x(1,:)); title('x_1[n]');
% subplot(2,2,2); plot(0:N, x(2,:)); title('x_2[n]');
% subplot(2,2,3); plot(0:N, x(3,:)); title('x_3[n]');
% subplot(2,2,4); plot(0:N, h);      title('h[n]');
% 
% h_filter = filter(b_tf, a_tf, [1; zeros(N,1)]);
% figure;
% plot(0:N, h, 'b', 0:N, h_filter, 'r--');
% legend('Estado', 'Filter');
% title('Comparaci贸n h[n]');
% 
% y1 = filter(b, a, v1);
% y2 = filter(b, a, v2);
% y3 = filter(b, a, v3);
% 
% v = repmat(v1, 25, 1);  % donde vp es una secuencia de longitud 4
% y = filter(b, a, v);

%-------------------Ejercicio 4---------------------
% 
% % N = 67;
% % fs = 8000;
% % Wn = 1000 / (fs/2); 
% % b = fir1(N, Wn, rectwin(N+1));
% % fvtool(b, 1); % fdatool para visualizacion
% 
% 
% % N = 67;
% % fs = 8000;
% % fc = 3000; 
% % Wn = fc/(fs/2); 
% % b = fir1(N, Wn, 'high', rectwin(N+2));
% % fvtool(b, 1);
% 
% % N = 67;
% % fs = 8000;
% % Wn = [1000 2000] / (fs/2); 
% % b = fir1(N, Wn, 'bandpass', rectwin(N+1));
% % fvtool(b, 1);
% 
% % N = 67;
% % fs = 8000;
% % Wn = [1000 2000] / (fs/2); 
% % b = fir1(N, Wn, 'stop', rectwin(N+2));
% % fvtool(b, 1);


% N = 20;
% fs = 8000;
% Wn = 1000 / (fs/2); 
% b = fir1(N, Wn, rectwin(N+1));
% fvtool(b, 1); % fdatool para visualizacion
%

%---------------------------Ejercicio 5 ---------------------------




%---------------------------Ejercicio 7 ----------------------
% %Dise駉 de un filtro noch (elimina frecuencia muy fino
% %La frecuiencia de eliminacion se toma pi*fe/fs 
% %Por lo tanto theeta corte = 2/5 pi
% %El noch FIR se caracteriza por tener todos los polos en el origen y con
% %ceros conjugados
% % K=? El modulo del filtro tiene que valer 1 en continua
% % =>k=1/(2(1-cos(Theeta corte))) 
% %El dise駉 del filtro cuenta con buscar la constante y a su vez saber bien las propiedades 
% 
% fe=50; %hz
% fs=250; 
% Limites = fs/2; %(tomo la mitad del plano, porque ser韆 como mi limite)
% 
% theetacorte = pi*2/5; %Sale de la regla de 3 mencionada arriba
% h = [1/(2*(1-cos(theetacorte))) -cos(theetacorte)./(1-cos(theetacorte)) 1/(2*(1-cos(theetacorte))) ];
%     
% [H w] = freqz(h,1,512);
% 
% %Grafico en frecuencia
% 
% f=w*(fs/(2*pi));
% % Abrir el archivo del ecg
% fid = fopen('ecg_nt_Ejercicio_7.txt', 'r');
% if fid == -1
%     error('No se pudo abrir el archivo.');
% end
% 
% % Leer datos (asume una columna de datos num閞icos)
% ecg = fscanf(fid, '%f\n');
% 
% % Cerrar archivo
% fclose(fid);
% 
% t = (0:length(ecg)-1)/fs;
% 
% %Filtracion de la se馻l
% 
% y = filter(h,1,ecg);
% % Graficar se馻l ECG
% figure(1)
% subplot(2,1,1)
% plot(t, ecg);
% xlabel('Tiempo (s)');
% ylabel('Amplitud');
% title('Se馻l ECG desde archivo con fopen');
% grid on;
% subplot(2,1,2)
% plot(t,y)
% %Graficacion del modulo del filtro notch 
% figure(2)
% plot(f,abs(H))



%-------------------Ejercicio 8---------------------

% %-------------------Ejercicio 10---------------------
% 
% fs = 1*10^3; %hz
% wc=2*pi*60; %rad
% num = [wc]; 
% den = [1 wc];
% Hs=tf(num,den); 
% bode(Hs);
% [numz,denz]=bilinear(num,den,fs);
% [TFS w]=freqz(numz,denz,512);
% Graficador_Freqz(TFS,w,'Respuesta Filtro RC');
% 


%-------------------Ejercicio 11---------------------
% 
% fs = 1*10^3; %hz
% wc=2*pi*60; %rad
% num = [wc]; 
% den = [1 wc];
% Hs=tf(num,den); 
% % bode(Hs);
% [numz,denz]=bilinear(num,den,fs);
% [TFS w]=freqz(numz,denz,512);
% Graficador_Freqz(TFS,w,'Respuesta Filtro RC');
% 
% t=0:1/fs:(512-1)*(1/fs);
% xt = 2*sin(2*pi*50*t)+ sin(2*pi*300*t);
% figure(2)
% subplot(2,1,1)
% plot(t,xt)
% yt = filter(numz,denz,xt)
% subplot(2,1,2)
% plot(yt);

%-------------------Ejercicio 12---------------------
% fs = 200; %hz
% w0=2*pi*50; %(1/s)
% Q0=100;
% num = [1 0 w0^2]; 
% den = [1 w0/Q0 w0^2];
% Hs=tf(num,den); 
% %bode(Hs);
% [numz,denz]=bilinear(num,den,fs);
% [TFS w]=freqz(numz,denz,512);
% Graficador_Freqz(TFS,w,'Respuesta Filtro RC');

%-------------------Ejercicio 13 - 17 paso de filtros ---------------------
% 
 % Especificaciones
fs = 1000;           % Frecuencia de muestreo en Hz
fp = 100;            % Frecuencia de pasabanda en Hz
fsb = 200;           % Frecuencia de detenci髇 en Hz
Rp = 10;              % Atenuaci髇 en pasabanda (dB)
Rs = 100;            % Atenuaci髇 en banda de detenci髇 (dB)

% Normalizaci髇 de frecuencias (0 a 1, donde 1 corresponde a fs/2)
Wp = 1 *2*pi;    
Ws = 1/2 *2*pi;   

% Orden de Butterworth
[ButOrd, Wn_butter] = buttord(Wp, Ws, Rp, Rs,'s');

% Orden de Chebyshev tipo I
[ChevOrd, Wn_cheby] = cheb1ord(Wp, Ws, Rp, Rs,'s');

% Orden de El韕tico
[EliOrd, Wn_ellip] = ellipord(Wp, Ws, Rp, Rs,'s');

% Dise駉 de filtros con sus 髍denes respectivos
[bButter, aButter] = butter(ButOrd, Wn_butter,'s');
[bCheby, aCheby]   = cheby1(ChevOrd, Rp, Wn_cheby,'s');
[bElipt, aElipt]   = ellip(EliOrd, Rp, Rs, Wn_ellip,'s');

%Pasaaltoscheby
[HighbChevy,HighaCheby]=lp2hp(bCheby,aCheby,2*pi*1000);
[numz,denz]=bilinear(HighbChevy,HighaCheby,fs);
 
% Frecuencia para evaluar la respuesta
[Hbutter, f] = freqz(bButter, aButter, 1024, fs);
Hcheby = freqz(bCheby, aCheby, 1024, fs);
Hellip = freqz(bElipt, aElipt, 1024, fs);

%Cheby pasa altos
[HighHcheby,w] = freqz(numz, denz, 1024, fs);

plot(w,abs(HighHcheby))

% % Subplot de las respuestas en frecuencia
% figure;
% subplot(3,1,1);
% plot(f, 20*log10(abs(Hbutter)));
% grid on;
% title(['Filtro Butterworth (Orden = ' num2str(ButOrd) ')']);
% xlabel('Frecuencia (Hz)');
% ylabel('Magnitud (dB)');
% ylim([-120 5]);
% 
% subplot(3,1,2);
% plot(f, 20*log10(abs(Hcheby)));
% grid on;
% title(['Filtro Chebyshev I (Orden = ' num2str(ChevOrd) ')']);
% xlabel('Frecuencia (Hz)');
% ylabel('Magnitud (dB)');
% ylim([-120 5]);
% 
% subplot(3,1,3);
% plot(f, 20*log10(abs(Hellip)));
% grid on;
% title(['Filtro El韕tico (Orden = ' num2str(EliOrd) ')']);
% xlabel('Frecuencia (Hz)');
% ylabel('Magnitud (dB)');
% ylim([-120 5]);

