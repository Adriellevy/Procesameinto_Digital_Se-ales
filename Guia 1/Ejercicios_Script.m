%Ejercicio 1  
 
% n = -10:10; % Rango de valores para n
% delta = (n == 0); % Genera el impulso de Kronecker
% 
% subplot(3,1,1)
% stem(n, delta, 'filled'); % Representación gráfica discreta
% xlabel('n');
% ylabel('\delta[n]');
% title('Impulso de Kronecker');
% grid on;
% 
% subplot(3,1,2)
% [x,y] = Escalon(10,1);
% stem(x,y,'filled'); % Representación gráfica discreta
% xlabel('x');
% ylabel('y[x]');
% title('Escalon');
% grid on;
% 
% 
% subplot(3,1,3)
% [x,y] = Rampa(10,1,1);
% stem(x,y,'filled'); % Representación gráfica discreta
% xlabel('x');
% ylabel('y[x]');
% title('Escalon');
% grid on;

%Ejercicio 2
 
% subplot(2,1,1)
% [x,y] = Cierra(10,1,2,4);
% stem(x,y,'filled'); % Representación gráfica discreta
% xlabel('x');
% ylabel('y[x]');
% title('Cierra');
% grid on;
% 
% subplot(2,1,2)
% [x,y] = Onda_Cuadrada(10,1,2,4,50);
% stem(x, y, 'filled');
% xlabel('Tiempo (s)');
% ylabel('Amplitud');
% grid on;

%Ejercicio 3
% 
% figure(1)
% [x,func]=Funcion_senoidal_w(-15,25,3,pi*1/10,0);
% stem(x, func, 'filled');
% xlabel('Tiempo (s)');
% ylabel('Amplitud');
% grid on;



%Ejercicio 4
clc; clear;

% % Parámetros de la señal

% L = 10; 
% A = 1;         % Amplitud
% phi = 0;       % Fase en radianes
% f0 = 5;        % Frecuencia de la señal en Hz
% fs = f0/3;
% % === Definición del tiempo ===
% T0 = 1 / f0;   % Período de la señal
% t_cont = 0:T0/100:L*T0;  % Tiempo continuo (para trazar bien la senoidal)
% t_samp = 0:1/fs:L*T0;    % Tiempo de muestreo
% 
% % === Generación de señales ===
% x_cont = A * sin(2 * pi * f0 * t_cont + phi);  % Señal continua
% x_samp = A * sin(2 * pi * f0 * t_samp + phi);  % Señal muestreada
% 
% % === Gráficos ===
% figure;
% plot(t_cont, x_cont, 'b', 'LineWidth', 1.5);  % Señal continua
% hold on;
% stem(t_samp, x_samp, 'ro', 'MarkerSize', 4, 'LineWidth', 1.5);  % Señal muestreada
% hold off;
% 
% % Configuración del gráfico
% xlabel('Tiempo (s)');
% ylabel('Amplitud');
% title(['Señal Senoidal y su Versión Muestreada (fs = ', num2str(fs), ' Hz)']);
% grid on;
% legend('Señal Continua', 'Muestreada');

%Ejercicio 5
%Se hace uso de la funcion impz en mi Funcion_Respuesta_discreta

A = [1 0.7 -0.45 -0.6]; %terminos que acompañan a las y
B = [0.8 -0.44 0.16 0.02];  %terminos que acompañan a las x
N = 50;
Funcion_Respuesta_discreta(A,B,N)


%Ejercicio 6
%Se hace uso de la funcion impz
%Inciso a

%Inciso b
A = [1 0 0.9]; %terminos que acompañan a las y
B = [.3 .6 .3];  %terminos que acompañan a las x
N = 163;
figure(1)
Funcion_Respuesta_discreta(A,B,N)

A = [1 1.8*cos(pi/16) 0.81]; 
B = [1 1/2 0];
figure(2)
Funcion_Respuesta_discreta(A,B,N)