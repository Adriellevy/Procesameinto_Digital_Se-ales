%-------------------Ejercicio 1---------------------
% a = [1, -0.8741, 0.9217, -0.2672];
% b = [0.1866, 0.2036, 0.2036, 0.1866];
% 
% % Representación en espacio de estados
% [A, B, C, D] = tf2ss(b, a);
% 
% % Expansión en fracciones simples
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
% % Comparación con fracciones simples
% h_analitica = zeros(100, 1);
% for i = 1:length(r)
%     h_analitica = h_analitica + real(r(i) * (p(i).^(0:99)'));
% end
% 
% % Comparación visual
% figure;
% plot(0:99, h, 'b', 0:99, h_analitica, 'r--');
% legend('Respuesta filter', 'Respuesta fracciones simples');
% title('Comparación respuesta al impulso');

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
% autovalores = eig(A);  % deberían coincidir
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
% title('Comparación h[n]');
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
% % fvtool(b, 1); % fdatool para visualización
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
% 
% % 
% % N = 20;
% % fs = 8000;
% % Wn = 1000 / (fs/2); 
% % b = fir1(N, Wn, rectwin(N+1));
% % fvtool(b, 1); % fdatool para visualización