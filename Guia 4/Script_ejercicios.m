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

%plot(mse)

%Automatico
k=0:0.1:10;
N = length(k); % cantidad de muestras
x = sin(2*pi*k);       % entrada
d = 2*cos(2*pi*k); 
w0 = zeros(1, N);       
w1 = zeros(1, N);
y = zeros(1, N);
e = zeros(1, N);  
mu = 0.5;
for k = 2:N
    y(k) = w0(k-1)*x(k) + w1(k-1)*x(k-1);    % salida del filtro
    e(k) = d(k) - y(k);                      % error
    w0(k) = w0(k-1) + 2*mu*e(k)*x(k);        % actualizar peso 0
    w1(k) = w1(k-1) + 2*mu*e(k)*x(k-1);      % actualizar peso 1
end

figure(1)
subplot(2,1,1)
plot(w0); hold on; plot(w1); legend('w0','w1');
title('Evolución de los pesos');
hold off;
subplot(2,1,2)
plot(d); hold on; plot(y); legend('d','y');
title('Comparacion funciones');
hold off;